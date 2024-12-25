#include "pipeline.hpp"
#include "IShader.hpp"
#include "camera.hpp"
#include "framebuffer.hpp"
#include <cassert>
#include <algorithm>
#include <cmath>

namespace AR {
	// Constructor
	Pipeline::Pipeline()
		: m_Shader(nullptr)
		, m_Camera(nullptr)
		, m_Framebuffer(nullptr)
	{
	}

	// Setters
	void Pipeline::setShader(IShader* shader) {
		m_Shader = shader;
	}

	void Pipeline::setCamera(const Camera* cam) {
		m_Camera = cam;
	}

	void Pipeline::setFramebuffer(Framebuffer* fb) {
		m_Framebuffer = fb;
	}

	// Get viewport matrix
	mat4f Pipeline::getViewportMat()
	{
		return m_Camera->getViewportMatrix();
	}

	// Draw a mesh
	void Pipeline::drawMesh(const mat4f& modelMatrix, const Mesh& mesh) {
		// Check if everything is set up
		if (!m_Shader || !m_Camera || !m_Framebuffer) return;

		// Precompute transformation matrices
		mat4f view = m_Camera->getViewMatrix();
		mat4f proj = m_Camera->getProjectionMatrix();
		mat4f viewProj = proj * view;
		mat4f mvp = viewProj * modelMatrix;

		// Set shader uniforms
		m_Shader->model = modelMatrix;
		m_Shader->viewProj = viewProj;
		m_Shader->mvp = mvp;
		m_Shader->viewportMat = m_Camera->getViewportMatrix();
		m_Shader->cameraPosition = m_Camera->getPosition();

		// Get mesh data
		const Vertex* vertices = mesh.getVertices().data();
		const std::vector<Face>& faces = mesh.getFaces();
		const auto& groups = mesh.getMaterialGroups();

		// Iterate over material groups
		for (const auto& group : groups) {
			// Set material
			const Material* material = mesh.getMaterial(group.materialName);
			m_Shader->material = material;

			// Iterate over faces in the group
			for (size_t i = group.startIndex; i < group.startIndex + group.faceCount; ++i) {
				const auto& face = faces[i];
				assert(face.vertexIndices.size() == 3);

				// Get vertex data
				const Vertex& v0Data = vertices[face.vertexIndices[0]];
				const Vertex& v1Data = vertices[face.vertexIndices[1]];
				const Vertex& v2Data = vertices[face.vertexIndices[2]];

				// Vertex shader stage
				Vec4f clipCoords[3];
				clipCoords[0] = mvp * toVec4f(v0Data.position);
				clipCoords[1] = mvp * toVec4f(v1Data.position);
				clipCoords[2] = mvp * toVec4f(v2Data.position);

				std::vector<std::array<std::pair<Vertex, Vec4f>, 3>> clippedTriangles = clipTriangle({
					std::make_pair<>(v0Data, clipCoords[0]),
					std::make_pair<>(v1Data, clipCoords[1]),
					std::make_pair<>(v2Data, clipCoords[2]),
					});

				for (auto& tri : clippedTriangles)
				{
					clipCoords[0] = m_Shader->vertex(tri[0].first, 0);
					clipCoords[1] = m_Shader->vertex(tri[1].first, 1);
					clipCoords[2] = m_Shader->vertex(tri[2].first, 2);

					rasterizeTriangle(clipCoords);
				}
			}
		}
	}

	//-----------------------------------------------
	//    Clipping Implementation
	//-----------------------------------------------

	// Clip a triangle against the six frustum planes
	std::vector<std::array<std::pair<Vertex, Vec4f>, 3>> Pipeline::clipTriangle(const std::array<std::pair<Vertex, Vec4f>, 3>& tri) {
		std::vector<std::pair<Vertex, Vec4f>> polygon(tri.begin(), tri.end());

		// Clip against each plane
		for (int planeIndex = 0; planeIndex < 6; ++planeIndex) {
			polygon = clipAgainstPlane(polygon, planeIndex);
			if (polygon.size() < 3) {
				break; // Triangle is fully clipped
			}
		}

		// Triangulate the clipped polygon (if necessary)
		std::vector<std::array<std::pair<Vertex, Vec4f>, 3>> outTriangles;
		if (polygon.size() >= 3) {
			for (size_t i = 1; i < polygon.size() - 1; ++i) {
				outTriangles.push_back({ polygon[0], polygon[i], polygon[i + 1] });
			}
		}

		return outTriangles;
	}

	bool Pipeline::insidePlane(const Vec4f& v, int plane) {
		switch (plane) {
		case 0: return (v.x + v.w) >= 0; // x >= -w
		case 1: return (v.w - v.x) >= 0; // x <=  w
		case 2: return (v.y + v.w) >= 0; // y >= -w
		case 3: return (v.w - v.y) >= 0; // y <=  w
		case 4: return (v.z + v.w) >= 0; // z >= -w
		case 5: return (v.w - v.z) >= 0; // z <=  w
		default: return false;
		}
	}

	std::pair<Vertex, Vec4f> Pipeline::interpolateVertices(std::pair<Vertex, Vec4f> v0, std::pair<Vertex, Vec4f> v1, float t_Point)
	{
		std::pair<Vertex, Vec4f> out;
		out.first.position = v0.first.position + (v1.first.position - v0.first.position) * t_Point;
		out.first.normal = v0.first.normal + (v1.first.normal - v0.first.normal) * t_Point;
		out.first.uv = v0.first.uv + (v1.first.uv - v0.first.uv) * t_Point;
		out.first.tangent = v0.first.tangent + (v1.first.tangent - v0.first.tangent) * t_Point;
		out.first.bitangent = v0.first.bitangent + (v1.first.bitangent - v0.first.bitangent) * t_Point;

		out.second = v0.second + (v1.second - v0.second) * t_Point;

		return out;
	}

	float Pipeline::intersectPlane(const Vec4f& v1, const Vec4f& v2, int plane) {
		// Helper lambda for plane distance calculation
		auto distFunc = [&](const Vec4f& v, int pl) {
			switch (pl) {
			case 0: return  v.x + v.w;
			case 1: return  v.w - v.x;
			case 2: return  v.y + v.w;
			case 3: return  v.w - v.y;
			case 4: return  v.z + v.w;
			case 5: return  v.w - v.z;
			default: return 0.0f;
			}
			};

		float d1 = distFunc(v1, plane);
		float d2 = distFunc(v2, plane);

		// Avoid division by zero
		float denom = (d1 - d2);
		if (std::fabs(denom) < 1e-7f) {
			return 0.5f; // Fallback: midpoint
		}
		float t = d1 / denom;
		return t;
	}

	// Clip a polygon against a single plane
	std::vector<std::pair<Vertex, Vec4f>> Pipeline::clipAgainstPlane(const std::vector<std::pair<Vertex, Vec4f>>& poly, int plane) {
		std::vector<std::pair<Vertex, Vec4f>> out;
		if (poly.empty()) return out;

		// Sutherland-Hodgman algorithm
		for (size_t i = 0; i < poly.size(); ++i) {
			const std::pair<Vertex, Vec4f>& currentPair = poly[i];
			const std::pair<Vertex, Vec4f>& nextPair = poly[(i + 1) % poly.size()];

			const Vec4f& currentPos = currentPair.second;
			const Vec4f& nextPos = nextPair.second;

			bool cInside = insidePlane(currentPos, plane);
			bool nInside = insidePlane(nextPos, plane);

			if (cInside && nInside) {
				out.push_back(nextPair); // Both inside
			}
			else if (cInside && !nInside) {
				float t = intersectPlane(currentPos, nextPos, plane);
				out.push_back(interpolateVertices(currentPair, nextPair, t)); // Current inside, next outside
			}
			else if (!cInside && nInside) {
				float t = intersectPlane(currentPos, nextPos, plane);
				out.push_back(interpolateVertices(currentPair, nextPair, t)); // Current inside, next outside
				out.push_back(nextPair);
			}
		}

		return out;
	}

	//-----------------------------------------------
	//    Rasterization
	//-----------------------------------------------

	// Calculate barycentric coordinates
	Vec3f Pipeline::barycentric(const Vec2f& A, const Vec2f& B, const Vec2f& C, const Vec2f& P) const {
		float areaABC = (B.y - C.y) * (A.x - C.x) + (C.x - B.x) * (A.y - C.y); // Area of the triangle ABC
		float areaPBC = (B.y - C.y) * (P.x - C.x) + (C.x - B.x) * (P.y - C.y); // Area of the triangle PBC
		float areaPCA = (C.y - A.y) * (P.x - C.x) + (A.x - C.x) * (P.y - C.y); // Area of the triangle PCA

		float alpha = areaPBC / areaABC; // Alpha (weight for A)
		float beta = areaPCA / areaABC;  // Beta (weight for B)
		float gamma = 1.0f - alpha - beta; // Gamma (weight for C)

		if (std::abs(areaABC) < 1e-6f) {
			return Vec3f(-1.0f, 1.0f, 1.0f); // Degenerate triangle
		}

		return Vec3f(alpha, beta, gamma);
	}

	// Rasterize a triangle
	void Pipeline::rasterizeTriangle(const Vec4f clip[3]) {
		int width = m_Framebuffer->getWidth();
		int height = m_Framebuffer->getHeight();

		// Perspective division
		Vec3f ndc[3];
		for (int i = 0; i < 3; ++i) {
			float w = clip[i].w;
			// Prevent division by zero.
			if (std::abs(w) < 1e-6f)
				w = std::copysign(1e-6f, w);
			ndc[i] = Vec3f(clip[i].x / w, clip[i].y / w, clip[i].z / w);
		}

		// NDC to screen space
		Vec2i screen[3];
		for (int i = 0; i < 3; ++i) {
			screen[i].x = static_cast<int>((ndc[i].x + 1.0f) * 0.5f * width);
			screen[i].y = static_cast<int>((ndc[i].y + 1.0f) * 0.5f * height);
		}

		// Bounding box
		Vec2i bboxMin(width - 1, height - 1);
		Vec2i bboxMax(0, 0);
		for (int i = 0; i < 3; ++i) {
			bboxMin.x = std::max(0, std::min(bboxMin.x, screen[i].x));
			bboxMin.y = std::max(0, std::min(bboxMin.y, screen[i].y));
			bboxMax.x = std::min(width - 1, std::max(bboxMax.x, screen[i].x));
			bboxMax.y = std::min(height - 1, std::max(bboxMax.y, screen[i].y));
		}

		// Rasterization loop
		for (int y = bboxMin.y; y <= bboxMax.y; ++y) {
			for (int x = bboxMin.x; x <= bboxMax.x; ++x) {
				Vec2f P(static_cast<float>(x), static_cast<float>(y));
				Vec3f bcScreen = barycentric(Vec2f(screen[0].x, screen[0].y), Vec2f(screen[1].x, screen[1].y), Vec2f(screen[2].x, screen[2].y), P);

				// Check if the point is inside the triangle
				if (bcScreen.x < 0 || bcScreen.y < 0 || bcScreen.z < 0) continue;

				// Interpolate depth
				float z = ndc[0].z * bcScreen.x + ndc[1].z * bcScreen.y + ndc[2].z * bcScreen.z;

				// Depth test
				if (z < m_Framebuffer->getDepth(x, y)) {
					// Fragment shader
					Vec4f colorOut;
					if (!m_Shader->fragment(bcScreen, colorOut)) {
						// Set pixel color and depth
						m_Framebuffer->setDepth(x, y, z);
						m_Framebuffer->setPixel(x, y, {
							static_cast<uint8_t>(std::clamp(colorOut.x, 0.0f, 1.0f) * 255),
							static_cast<uint8_t>(std::clamp(colorOut.y, 0.0f, 1.0f) * 255),
							static_cast<uint8_t>(std::clamp(colorOut.z, 0.0f, 1.0f) * 255),
							static_cast<uint8_t>(std::clamp(colorOut.w, 0.0f, 1.0f) * 255)
							});
					}
				}
			}
		}
	}
}