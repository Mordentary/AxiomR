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
#include "tracy/Tracy.hpp"

namespace AR {
	Pipeline::Pipeline()
		: m_Shader(nullptr)
		, m_Camera(nullptr)
		, m_Framebuffer(nullptr)
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

			std::vector<std::vector<std::pair<Vertex, glm::vec4>>> allClippedVertices;
			std::vector<std::vector<std::array<std::pair<Vertex, glm::vec4>, 3>>> allClippedTriangles;

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

					std::vector<std::vector<std::pair<Vertex, glm::vec4>>> threadClippedVertices;
					std::vector<std::vector<std::array<std::pair<Vertex, glm::vec4>, 3>>> threadClippedTriangles;
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

						std::vector<std::pair<Vertex, glm::vec4>> clippedVertices;
						clippedVertices.reserve(15); // max possible after clipping

						clippedVertices.emplace_back(v0Data, clipCoords[0]);
						clippedVertices.emplace_back(v1Data, clipCoords[1]);
						clippedVertices.emplace_back(v2Data, clipCoords[2]);

						clipTriangle(clippedVertices);

						std::vector<std::array<std::pair<Vertex, glm::vec4>, 3>> clippedTriangles;
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
						std::lock_guard<std::mutex> lock(m_Mutex);
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

			// Rasterize all clipped triangles
			for (auto& clippedTriangles : allClippedTriangles) {
				for (auto& tri : clippedTriangles) {
					glm::vec4 clipCoords[3];
					clipCoords[0] = m_Shader->vertex(tri[0].first, 0);
					clipCoords[1] = m_Shader->vertex(tri[1].first, 1);
					clipCoords[2] = m_Shader->vertex(tri[2].first, 2);

					rasterizeTriangle(clipCoords);
				}
			}
		}
	}

	//-------------------------------------------------------------
	// Clipping Implementation
	//-------------------------------------------------------------
	void Pipeline::clipTriangle(std::vector<std::pair<Vertex, glm::vec4>>& clippedVertices) {
		ZoneScoped;

		// Clip against each plane
		for (int planeIndex = 0; planeIndex < 6; ++planeIndex) {
			clipAgainstPlane(clippedVertices, planeIndex);
			if (clippedVertices.size() < 3) {
				break; // Triangle is fully clipped
			}
		}
	}

	bool Pipeline::insidePlane(const glm::vec4& v, int plane) {
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

	std::pair<Vertex, glm::vec4> Pipeline::interpolateVertices(
		std::pair<Vertex, glm::vec4> v0,
		std::pair<Vertex, glm::vec4> v1,
		float t_Point)
	{
		ZoneScoped;
		std::pair<Vertex, glm::vec4> out;

		// Interpolate vertex attributes
		out.first.position = v0.first.position + (v1.first.position - v0.first.position) * t_Point;
		out.first.normal = v0.first.normal + (v1.first.normal - v0.first.normal) * t_Point;
		out.first.uv = v0.first.uv + (v1.first.uv - v0.first.uv) * t_Point;
		out.first.tangent = v0.first.tangent + (v1.first.tangent - v0.first.tangent) * t_Point;
		out.first.bitangent = v0.first.bitangent + (v1.first.bitangent - v0.first.bitangent) * t_Point;

		// Interpolate clip space position
		out.second = v0.second + (v1.second - v0.second) * t_Point;

		return out;
	}

	float Pipeline::intersectPlane(const glm::vec4& v1, const glm::vec4& v2, int plane) {
		// Helper lambda for plane distance calculation
		auto distFunc = [&](const glm::vec4& v, int pl) {
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

	// Clip a polygon against a single plane, store the result in 'poly'
	void Pipeline::clipAgainstPlane(std::vector<std::pair<Vertex, glm::vec4>>& poly, int plane) {
		ZoneScoped;
		if (poly.empty()) return;

		std::vector<std::pair<Vertex, glm::vec4>> tempOut;
		tempOut.reserve(poly.size() + 3);

		// Sutherland-Hodgman algorithm
		for (size_t i = 0; i < poly.size(); ++i) {
			const std::pair<Vertex, glm::vec4>& currentPair = poly[i];
			const std::pair<Vertex, glm::vec4>& nextPair = poly[(i + 1) % poly.size()];

			const glm::vec4& currentPos = currentPair.second;
			const glm::vec4& nextPos = nextPair.second;

			bool cInside = insidePlane(currentPos, plane);
			bool nInside = insidePlane(nextPos, plane);

			if (cInside && nInside) {
				// Both inside
				tempOut.push_back(nextPair);
			}
			else if (cInside && !nInside) {
				// Current inside, next outside
				float t = intersectPlane(currentPos, nextPos, plane);
				tempOut.push_back(interpolateVertices(currentPair, nextPair, t));
			}
			else if (!cInside && nInside) {
				// Current outside, next inside
				float t = intersectPlane(currentPos, nextPos, plane);
				tempOut.push_back(interpolateVertices(currentPair, nextPair, t));
				tempOut.push_back(nextPair);
			}
			// else both outside -> discard
		}

		poly.swap(tempOut);
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
		ZoneScoped;
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
					if (!m_Shader->fragment(bcScreen, colorOut)) {
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