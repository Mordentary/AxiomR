#include "pipeline.hpp"
#include "IShader.hpp"
#include "camera.hpp"
#include "framebuffer.hpp"
#include <cassert>
#include <algorithm>
#include <cmath>
#include <thread>
#include <vector>
#include <mutex>
#include <array>
#include "tracy/Tracy.hpp"
#include <tiled_pipeline.hpp>

namespace AR {
	Pipeline::Pipeline(Camera* cam, Framebuffer* fb) :
		m_Camera(cam),
		m_Framebuffer(fb)
	{
	}

	void Pipeline::setShader(IShader* shader) {
		m_Shader = shader;
	}

	void Pipeline::setCamera(Camera* cam) {
		m_Camera = cam;
	}

	void Pipeline::setFramebuffer(Framebuffer* fb) {
		m_Framebuffer = fb;
	}

	glm::mat4 Pipeline::getViewportMat()
	{
		return m_Camera->getViewportMatrix();
	}

	void Pipeline::drawMesh(const glm::mat4& modelMatrix, const Mesh& mesh) {
		if (!m_Shader || !m_Camera || !m_Framebuffer) return;

		// Precompute transformation matrices
		const glm::mat4& view = m_Camera->getViewMatrix();
		const glm::mat4& proj = m_Camera->getProjectionMatrix();
		const glm::mat4& viewProj = m_Camera->getViewProjectionMatrix();
		glm::mat4 mvp = viewProj * modelMatrix;

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

		for (const auto& group : groups) {
			const Material* material = mesh.getMaterial(group.materialName);
			m_Shader->material = material;

			std::vector<std::vector<ClippedVertex>> allClippedVertices;
			std::vector<std::vector<std::array<ClippedVertex, 3>>> allClippedTriangles;

			std::vector<std::thread> threads;
			int numThreads = std::thread::hardware_concurrency();
			if (numThreads == 0) {
				numThreads = 1;
			}

			int facesPerThread = (group.faceCount + numThreads - 1) / numThreads;

			for (int t = 0; t < numThreads; ++t) {
				threads.emplace_back([&, t]() {
					int start = static_cast<int>(group.startIndex) + t * facesPerThread;
					int end = std::min(start + facesPerThread,
						static_cast<int>(group.startIndex + group.faceCount));

					std::vector<std::vector<ClippedVertex>> threadClippedVertices;
					std::vector<std::vector<std::array<ClippedVertex, 3>>> threadClippedTriangles;
					threadClippedVertices.reserve(end - start);
					threadClippedTriangles.reserve(end - start);

					for (int i = start; i < end; ++i) {
						const auto& face = faces[i];
						assert(face.vertexIndices.size() == 3);

						const Vertex& v0Data = vertices[face.vertexIndices[0]];
						const Vertex& v1Data = vertices[face.vertexIndices[1]];
						const Vertex& v2Data = vertices[face.vertexIndices[2]];

						glm::vec4 clipCoords[3];
						clipCoords[0] = mvp * glm::vec4(v0Data.position, 1.0f);
						clipCoords[1] = mvp * glm::vec4(v1Data.position, 1.0f);
						clipCoords[2] = mvp * glm::vec4(v2Data.position, 1.0f);

						std::vector<ClippedVertex> clippedVertices;
						clippedVertices.reserve(15); // max possible after clipping

						clippedVertices.emplace_back(v0Data, clipCoords[0]);
						clippedVertices.emplace_back(v1Data, clipCoords[1]);
						clippedVertices.emplace_back(v2Data, clipCoords[2]);

						//clipTriangle(clippedVertices);

						std::vector<std::array<ClippedVertex, 3>> clippedTriangles;
						if (clippedVertices.size() >= 3) {
							clippedTriangles.reserve(clippedVertices.size() - 2);
							for (size_t j = 1; j < clippedVertices.size() - 1; ++j) {
								clippedTriangles.push_back({
									clippedVertices[0],
									clippedVertices[j],
									clippedVertices[j + 1]
									});
							}
						}

						threadClippedVertices.push_back(clippedVertices);
						threadClippedTriangles.push_back(clippedTriangles);
					}

					{
						//std::lock_guard<std::mutex> lock(m_Mutex);
						allClippedVertices.insert(allClippedVertices.end(),
							threadClippedVertices.begin(),
							threadClippedVertices.end());
						allClippedTriangles.insert(allClippedTriangles.end(),
							threadClippedTriangles.begin(),
							threadClippedTriangles.end());
					}
					});
			}

			for (auto& thread : threads) {
				thread.join();
			}

			//// Rasterize all clipped triangles
			//for (auto& clippedTriangles : allClippedTriangles) {
			//    for (auto& tri : clippedTriangles) {
			//        glm::vec4 clipCoords[3];
			//        clipCoords[0] = m_Shader->vertex(tri[0].vertex, 0);
			//        clipCoords[1] = m_Shader->vertex(tri[1].vertex, 1);
			//        clipCoords[2] = m_Shader->vertex(tri[2].vertex, 2);

			//        rasterizeTriangle(clipCoords);
			//    }
			//}
		}
	}

	glm::vec4 planeNormalForIndex(int planeIndex) {
		switch (planeIndex) {
		case 0: // left:  v.x + v.w >= 0
			return glm::vec4(1.0f, 0.0f, 0.0f, 1.0f);
		case 1: // right: v.w - v.x >= 0
			return glm::vec4(-1.0f, 0.0f, 0.0f, 1.0f);
		case 2: // bottom: v.y + v.w >= 0
			return glm::vec4(0.0f, 1.0f, 0.0f, 1.0f);
		case 3: // top: v.w - v.y >= 0
			return glm::vec4(0.0f, -1.0f, 0.0f, 1.0f);
		case 4: // near: v.z + v.w >= 0
			return glm::vec4(0.0f, 0.0f, 1.0f, 1.0f);
		case 5: // far:  v.w - v.z >= 0
			return glm::vec4(0.0f, 0.0f, -1.0f, 1.0f);
		default:
			return glm::vec4(0.0f); // Return a zero vector for invalid indices.
		}
	}
	//-------------------------------------------------------------
	// Clipping Implementation
	//-------------------------------------------------------------

	void Pipeline::clipTriangle(size_t& index,
		std::array<ClippedVertex, MAX_CLIPPED_VERTS>& outTris)
	{
		// Temporary storage for the next stage
		std::array<ClippedVertex, MAX_CLIPPED_VERTS> nextTris{};
		int nextCount = 0;

		// Clip against 6 planes
		for (int planeIndex = 0; planeIndex < 6; ++planeIndex) {
			nextCount = 0;  // reset for this pass

			// Process existing triangles in groups of 3
			for (int i = 0; i < index; i += 3) {
				ClippedVertex& v0 = outTris[i + 0];
				ClippedVertex& v1 = outTris[i + 1];
				ClippedVertex& v2 = outTris[i + 2];
				ClippedVertex v3;
				// up to 4 clipped vertices
				int result = clipTriangleSinglePlane(planeIndex, v0, v1, v2, v3);

				if (result == 3) {
					// just 1 triangle => 3 vertices
					if (nextCount + 3 <= MAX_CLIPPED_VERTS) {
						// Use normal copy, not std::move
						nextTris[nextCount + 0] = v0;
						nextTris[nextCount + 1] = v1;
						nextTris[nextCount + 2] = v2;
						nextCount += 3;
					}
				}
				else if (result == 4) {
					// => 2 triangles => 6 vertices
					if (nextCount + 6 <= MAX_CLIPPED_VERTS) {
						nextTris[nextCount + 0] = v0;
						nextTris[nextCount + 1] = v1;
						nextTris[nextCount + 2] = v2;

						nextTris[nextCount + 3] = v0; // same v0
						nextTris[nextCount + 4] = v2; // same v2
						nextTris[nextCount + 5] = v3;
						nextCount += 6;
					}
				}
				// else (0 => fully clipped), skip
			}

			index = nextCount;

			outTris.swap(nextTris);
			// Early exit if nothing is left
			if (index == 0) break;
		}
	}

	inline bool Pipeline::insidePlane(const glm::vec4& v, int plane) {
		ZoneScoped;
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

	inline ClippedVertex Pipeline::interpolateVertices(
		const ClippedVertex& v0,
		const ClippedVertex& v1,
		float t_Point)
	{
		ZoneScoped;
		ClippedVertex out;

		// Interpolate vertex attributes using glm::mix
		out.vertex.position = glm::mix(v0.vertex.position, v1.vertex.position, t_Point);
		out.vertex.normal = glm::normalize(glm::mix(v0.vertex.normal, v1.vertex.normal, t_Point));
		out.vertex.tangent = glm::normalize(glm::mix(v0.vertex.tangent, v1.vertex.tangent, t_Point));
		out.vertex.bitangent = glm::normalize(glm::mix(v0.vertex.bitangent, v1.vertex.bitangent, t_Point));

		out.clipPos = glm::mix(v0.clipPos, v1.clipPos, t_Point);

		// Perspective-correct interpolation for UVs
		float w0 = v0.clipPos.w;
		float w1 = v1.clipPos.w;

		glm::vec2 uv0_weighted = v0.vertex.uv * w0;
		glm::vec2 uv1_weighted = v1.vertex.uv * w1;
		glm::vec2 interpolated_weighted_uv = glm::mix(uv0_weighted, uv1_weighted, t_Point);

		float interpolated_w = glm::mix(w0, w1, t_Point);

		out.vertex.uv = interpolated_weighted_uv / interpolated_w;

		return out;
	}

	static float distFunc(const glm::vec4& v, int plane)
	{
		switch (plane) {
		case 0: return  v.x + v.w;
		case 1: return  v.w - v.x;
		case 2: return  v.y + v.w;
		case 3: return  v.w - v.y;
		case 4: return  v.z + v.w;
		case 5: return  v.w - v.z;
		default: return 0.0f;
		}
	};

	inline float Pipeline::intersectPlane(const glm::vec4& v1, const glm::vec4& v2, int plane) {
		ZoneScoped;
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

	/*casual-effects.com/research/McGuire2011Clipping/McGuire-Clipping.pdf*/
	int Pipeline::clipTriangleSinglePlane(
		int plane,
		ClippedVertex& v0,
		ClippedVertex& v1,
		ClippedVertex& v2,
		ClippedVertex& v3)
	{
		float d0 = distFunc(v0.clipPos, plane);
		float d1 = distFunc(v1.clipPos, plane);
		float d2 = distFunc(v2.clipPos, plane);

		// Check if all are outside (below)
		if (d0 < 0.f && d1 < 0.f && d2 < 0.f) {
			return 0; // fully clipped
		}

		// Check if all are inside (above)
		if (d0 >= 0.f && d1 >= 0.f && d2 >= 0.f) {
			v3 = v0; // dummy
			return 3;
		}

		auto inside = [&](float d) { return d >= 0.f; };

		// Rotate the triangle so that v0 is guaranteed inside
		//   if it's not, we cycle (v0->v1, v1->v2, v2->v0)
		if (inside(d1) && !inside(d0)) {
			std::swap(v0, v1);
			std::swap(d0, d1);
			std::swap(v1, v2);
			std::swap(d1, d2);
		}
		else if (inside(d2) && !inside(d1)) {
			// Rotate the other way
			std::swap(v2, v1);
			std::swap(d2, d1);
			std::swap(v1, v0);
			std::swap(d1, d0);
		}

		// Now v0 is definitely inside. Next we always find intersection along (v0->v2)
		float denom02 = (d0 - d2);
		float t02 = (std::fabs(denom02) < 1e-7f) ? 0.5f : (d0 / denom02);
		v3 = interpolateVertices(v0, v2, t02);

		// Determine if v1 is also inside
		bool v1Inside = (d1 >= 0.f);
		if (v1Inside) {
			// => 2 vertices inside => we become a quad
			// We also find intersection along (v1->v2)
			float denom12 = (d1 - d2);
			float t12 = (std::fabs(denom12) < 1e-7f) ? 0.5f : (d1 / denom12);
			v2 = interpolateVertices(v1, v2, t12);

			// Return 4 => (v0, v1, v2, v3)
			return 4;
		}
		else {
			// => Only 1 vertex inside => final shape is a smaller triangle
			// Clip the edge (v0->v1)
			float denom01 = (d0 - d1);
			float t01 = (std::fabs(denom01) < 1e-7f) ? 0.5f : (d0 / denom01);
			v1 = interpolateVertices(v0, v1, t01);
			v2 = v3;

			// v2 = v3 => we already found intersection on (v0->v2)
			return 3;
		}
	}

	// Clip a polygon against a single plane, store the result in 'poly'
	void Pipeline::clipAgainstPlane(std::pmr::vector<ClippedVertex>& poly, int plane, std::pmr::vector<ClippedVertex>& tempOut) {
		//ZoneScoped;
		//if (poly.empty()) return;

		//tempOut.clear();
		//tempOut.reserve(poly.size() + 3); // Reserve enough space

		//// Sutherland-Hodgman algorithm
		//for (size_t i = 0; i < poly.size(); ++i) {
		//	const ClippedVertex& currentV = poly[i];
		//	const ClippedVertex& nextV = poly[(i + 1) % poly.size()];

		//	const glm::vec4& currentPos = currentV.clipPos;
		//	const glm::vec4& nextPos = nextV.clipPos;

		//	bool currentInside = insidePlane(currentPos, plane);
		//	bool nextInside = insidePlane(nextPos, plane);

		//	float t;
		//	if (currentInside) {
		//		if (nextInside) {
		//			tempOut.emplace_back(nextV);
		//		}
		//		else {
		//			t = intersectPlane(currentPos, nextPos, plane);
		//			tempOut.emplace_back(interpolateVertices(currentV, nextV, t));
		//		}
		//	}
		//	else if (nextInside) {
		//		t = intersectPlane(currentPos, nextPos, plane);
		//		tempOut.emplace_back(interpolateVertices(currentV, nextV, t));
		//		tempOut.emplace_back(nextV);
		//	}
		//}

		//poly.swap(tempOut);
	}

	//-------------------------------------------------------------
	// Rasterization
	//-------------------------------------------------------------
	glm::vec3 Pipeline::barycentric(const glm::vec2& A, const glm::vec2& B,
		const glm::vec2& C, const glm::vec2& P) const
	{
		float areaABC = (B.y - C.y) * (A.x - C.x) + (C.x - B.x) * (A.y - C.y);
		float areaPBC = (B.y - C.y) * (P.x - C.x) + (C.x - B.x) * (P.y - C.y);
		float areaPCA = (C.y - A.y) * (P.x - C.x) + (A.x - C.x) * (P.y - C.y);

		// Avoid degenerate triangles
		if (std::abs(areaABC) < 1e-6f) {
			// Return an indicator that barycentric cannot be computed
			return glm::vec3(-1.0f, 1.0f, 1.0f);
		}

		float alpha = areaPBC / areaABC;
		float beta = areaPCA / areaABC;
		float gamma = 1.0f - alpha - beta;
		return glm::vec3(alpha, beta, gamma);
	}

	void Pipeline::rasterizeTriangle(const glm::vec4 clip[3]) {
		//ZoneScoped;
		int width = m_Framebuffer->getWidth();
		int height = m_Framebuffer->getHeight();

		// Perspective division to get NDC
		glm::vec3 ndc[3];
		for (int i = 0; i < 3; ++i) {
			float w = clip[i].w;
			if (std::abs(w) < 1e-6f) {
				w = std::copysign(1e-6f, w);
			}
			ndc[i] = glm::vec3(clip[i].x / w, clip[i].y / w, clip[i].z / w);
		}

		// Convert NDC to screen coordinates
		glm::ivec2 screen[3];
		for (int i = 0; i < 3; ++i) {
			screen[i].x = static_cast<int>((ndc[i].x + 1.0f) * 0.5f * width);
			screen[i].y = static_cast<int>((ndc[i].y + 1.0f) * 0.5f * height);
		}

		// Compute bounding box
		glm::ivec2 bboxMin(width - 1, height - 1);
		glm::ivec2 bboxMax(0, 0);
		for (int i = 0; i < 3; ++i) {
			bboxMin.x = std::max(0, std::min(bboxMin.x, screen[i].x));
			bboxMin.y = std::max(0, std::min(bboxMin.y, screen[i].y));
			bboxMax.x = std::min(width - 1, std::max(bboxMax.x, screen[i].x));
			bboxMax.y = std::min(height - 1, std::max(bboxMax.y, screen[i].y));
		}

		// Rasterize
		for (int y = bboxMin.y; y <= bboxMax.y; ++y) {
			for (int x = bboxMin.x; x <= bboxMax.x; ++x) {
				glm::vec2 P(static_cast<float>(x), static_cast<float>(y));

				glm::vec3 bcScreen = barycentric(
					glm::vec2(screen[0].x, screen[0].y),
					glm::vec2(screen[1].x, screen[1].y),
					glm::vec2(screen[2].x, screen[2].y),
					P
				);

				// If point is outside the triangle, skip
				if (bcScreen.x < 0.0f || bcScreen.y < 0.0f || bcScreen.z < 0.0f) {
					continue;
				}

				// Interpolate depth
				float z = ndc[0].z * bcScreen.x +
					ndc[1].z * bcScreen.y +
					ndc[2].z * bcScreen.z;

				// Depth test
				if (z < m_Framebuffer->getDepth(x, y)) {
					// Fragment shader
					glm::vec4 colorOut;
					// If fragment(...) returns false, we assume it wrote to colorOut
					//if (!m_Shader->fragment(bcScreen, colorOut))
					{
						m_Framebuffer->setDepth(x, y, z);
						m_Framebuffer->setPixel(x, y, {
							static_cast<uint8_t>(std::clamp(colorOut.r, 0.0f, 1.0f) * 255.0f),
							static_cast<uint8_t>(std::clamp(colorOut.g, 0.0f, 1.0f) * 255.0f),
							static_cast<uint8_t>(std::clamp(colorOut.b, 0.0f, 1.0f) * 255.0f),
							static_cast<uint8_t>(std::clamp(colorOut.a, 0.0f, 1.0f) * 255.0f)
							});
					}
				}
			}
		}
	}
}