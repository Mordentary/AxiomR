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
#include <atomic>
#include <immintrin.h>
#include <thread>
#include <queue>
#include <array>
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

			//// Rasterize all clipped triangles
			//for (auto& clippedTriangles : allClippedTriangles) {
			//    for (auto& tri : clippedTriangles) {
			//        glm::vec4 clipCoords[3];
			//        clipCoords[0] = m_Shader->vertex(tri[0].first, 0);
			//        clipCoords[1] = m_Shader->vertex(tri[1].first, 1);
			//        clipCoords[2] = m_Shader->vertex(tri[2].first, 2);

			//        rasterizeTriangle(clipCoords);
			//    }
			//}
		}
	}

	//-------------------------------------------------------------
	// Clipping Implementation
	//-------------------------------------------------------------
	void Pipeline::clipTriangle(std::vector<std::pair<Vertex, glm::vec4>>& clippedVertices) {
		//ZoneScoped;

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
		//ZoneScoped;
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
		//ZoneScoped;
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

	namespace
	{
		static inline float clampW(float w)
		{
			// Force a tiny nonzero magnitude if w is extremely close to zero
			const float tiny = 1e-6f;
			if (std::abs(w) < tiny) {
				w = (w < 0.0f) ? -tiny : tiny;
			}
			return w;
		}
	}
	// Thread-local buffers definition
	thread_local TiledPipeline::ThreadLocalBuffers TiledPipeline::t_buffers(TILE_SIZE);

	// ---------------- Triangle methods ----------------
	Triangle::Triangle(const glm::vec4* clipSpace, const Vertex* verts,
		int width, int height, const Material* mat)
		: material(mat)
	{
		for (int i = 0; i < 3; i++) {
			vertices[i] = verts[i];

			float w = clampW(clipSpace[i].w);
			float invW = 1.0f / w;

			// Z in NDC
			ndcZ[i] = clipSpace[i].z * invW;

			// Screen coordinates
			screenPos[i].x = ((clipSpace[i].x * invW) + 1.0f) * 0.5f * width;
			screenPos[i].y = ((clipSpace[i].y * invW) + 1.0f) * 0.5f * height;
		}
		// Triangle bounding box in screen space
		minX = std::min({ screenPos[0].x, screenPos[1].x, screenPos[2].x });
		minY = std::min({ screenPos[0].y, screenPos[1].y, screenPos[2].y });
		maxX = std::max({ screenPos[0].x, screenPos[1].x, screenPos[2].x });
		maxY = std::max({ screenPos[0].y, screenPos[1].y, screenPos[2].y });
	}

	Triangle::Triangle(const std::array<std::pair<Vertex, glm::vec4>, 3>& verts, int width, int height, const Material* mat)
		: material(mat)
	{
		for (size_t i = 0; i < verts.size(); ++i) {
			vertices[i] = verts[i].first;

			float w = clampW(verts[i].second.w);
			float invW = 1.0f / w;

			ndcZ[i] = verts[i].second.z * invW;

			screenPos[i].x = ((verts[i].second.x * invW) + 1.0f) * 0.5f * width;
			screenPos[i].y = ((verts[i].second.y * invW) + 1.0f) * 0.5f * height;
		}
		// Triangle bounding box in screen space
		minX = std::min({ screenPos[0].x, screenPos[1].x, screenPos[2].x });
		minY = std::min({ screenPos[0].y, screenPos[1].y, screenPos[2].y });
		maxX = std::max({ screenPos[0].x, screenPos[1].x, screenPos[2].x });
		maxY = std::max({ screenPos[0].y, screenPos[1].y, screenPos[2].y });
	}

	bool Triangle::isBackface() const {
		const glm::vec2& p0 = screenPos[0];
		const glm::vec2& p1 = screenPos[1];
		const glm::vec2& p2 = screenPos[2];

		float dx1 = p1.x - p0.x;
		float dy1 = p1.y - p0.y;
		float dx2 = p2.x - p0.x;
		float dy2 = p2.y - p0.y;

		float signedArea = dx1 * dy2 - dx2 * dy1;

		return signedArea < 0;
	}

	// ---------------- TiledPipeline methods ----------------

	TiledPipeline::TiledPipeline(int numThreads)
	{
		//ZoneScopedN("TiledPipeline Constructor");
		// Spawn threads
		for (int i = 0; i < numThreads; ++i) {
			m_Workers.emplace_back([this]() { workerThread(); });
		}
	}

	TiledPipeline::~TiledPipeline()
	{
		//ZoneScopedN("TiledPipeline Destructor");
		{
			std::unique_lock<std::mutex> lock(m_QueueMutex);
			m_ShouldExit = true;
		}
		m_WorkAvailable.notify_all();

		for (auto& thread : m_Workers) {
			if (thread.joinable())
				thread.join();
		}
	}

	void TiledPipeline::drawMesh(const glm::mat4& modelMatrix, const Mesh& mesh)
	{
		ZoneScopedN("TiledPipeline::drawMesh");

		if (!m_Shader || !m_Camera || !m_Framebuffer) return;

		const glm::mat4& viewProj = m_Camera->getViewProjectionMatrix();
		glm::mat4 mvp = viewProj * modelMatrix;

		m_Shader->model = modelMatrix;
		m_Shader->viewProj = viewProj;
		m_Shader->mvp = mvp;
		m_Shader->viewportMat = m_Camera->getViewportMatrix();
		m_Shader->cameraPosition = m_Camera->getPosition();

		// Prepare tiles and results

		{
			ZoneScopedN("Initialize Tiles");
			initializeTiles();
			m_TileResults.clear();
			m_TileResults.resize(m_Tiles.size(), TileResult(TILE_SIZE));
		}

		m_Triangles.clear();
		//m_Triangles.reserve(mesh.getFaces().size() * 2);

		const Vertex* vertices = mesh.getVertices().data();
		const std::vector<Face>& faces = mesh.getFaces();
		const auto& groups = mesh.getMaterialGroups();

		uint32_t framebufferWidth = m_Framebuffer->getWidth();
		uint32_t framebufferHeight = m_Framebuffer->getHeight();

		for (const auto& group : groups)
		{
			ZoneScopedN("Processing Material Group");
			const Material* material = mesh.getMaterial(group.materialName);
			m_Shader->material = material;

			// split clipping across concurrency to speed it up
			std::vector<std::thread> threads;
			int numThreads = std::thread::hardware_concurrency();
			if (numThreads == 0) { numThreads = 1; }

			int facesPerThread = (int)((group.faceCount + numThreads - 1) / numThreads);

			std::vector<std::vector<Triangle>> allClippedTriangles;
			allClippedTriangles.resize(numThreads);

			for (int t = 0; t < numThreads; ++t) {
				threads.emplace_back([&, t]() {
					ZoneScopedN("Clipping Thread");

					int start = (int)(group.startIndex) + t * facesPerThread;
					int end = std::min(start + facesPerThread,
						(int)(group.startIndex + group.faceCount));

					std::vector<Triangle> threadClippedTriangles;
					threadClippedTriangles.reserve((size_t)(end - start));

					for (int i = start; i < end; ++i) {
						const auto& face = faces[i];
						assert(face.vertexIndices.size() == 3);

						const Vertex& v0 = vertices[face.vertexIndices[0]];
						const Vertex& v1 = vertices[face.vertexIndices[1]];
						const Vertex& v2 = vertices[face.vertexIndices[2]];

						// Compute clip-space coords
						glm::vec4 clipCoords[3];
						clipCoords[0] = mvp * glm::vec4(v0.position, 1.0f);
						clipCoords[1] = mvp * glm::vec4(v1.position, 1.0f);
						clipCoords[2] = mvp * glm::vec4(v2.position, 1.0f);

						// Collect into a small vector for clipping
						std::vector<std::pair<Vertex, glm::vec4>> clippedVertices;
						clippedVertices.reserve(15);
						clippedVertices.emplace_back(v0, clipCoords[0]);
						clippedVertices.emplace_back(v1, clipCoords[1]);
						clippedVertices.emplace_back(v2, clipCoords[2]);

						// Frustum clipping in place
						clipTriangle(clippedVertices);

						// Triangulate any polygon produced by clipping
						if (clippedVertices.size() >= 3) {
							for (size_t j = 1; j < clippedVertices.size() - 1; ++j) {
								const auto& cv0 = clippedVertices[0];
								const auto& cv1 = clippedVertices[j];
								const auto& cv2 = clippedVertices[j + 1];
								Triangle newTri{ std::array{ cv0, cv1, cv2 },
										(int)framebufferWidth, (int)framebufferHeight,

										material };

								if (!newTri.isBackface()) {
									threadClippedTriangles.push_back(std::move(newTri));
								}
							}
						}
					}
					allClippedTriangles[t] = std::move(threadClippedTriangles);
					});
			}

			for (auto& th : threads) {
				th.join();
			}

			{
				ZoneScopedN("Flattening Clipped Triangles");
				size_t totalTriangles = 0;
				for (const auto& innerVec : allClippedTriangles) {
					totalTriangles += innerVec.size();
				}

				m_Triangles.reserve(totalTriangles);
				for (auto& innerVec : allClippedTriangles) {
					m_Triangles.insert(
						m_Triangles.end(),
						std::make_move_iterator(innerVec.begin()),
						std::make_move_iterator(innerVec.end())
					);
				}
			}

			{
				ZoneScopedN("Binning Triangles to Tiles");
				binTrianglesToTiles();
			}

			// Kick off tile processing
			m_NextTile = 0;
			m_CompletedTiles = 0;
			m_TotalTiles = m_Tiles.size();

			{
				std::unique_lock<std::mutex> lock(m_QueueMutex);
				m_HasWork = true;
			}
			m_WorkAvailable.notify_all();

			// Wait until all tiles processed
			{
				ZoneScopedN("Waiting for Tile Processing");
				std::unique_lock<std::mutex> lock(m_QueueMutex);
				m_WorkComplete.wait(lock, [this]() {
					return (m_CompletedTiles == m_TotalTiles);
					});
			}

			// Turn off the "work" signal
			m_HasWork = false;

			{
				// Merge tile results into the main framebuffer
				ZoneScopedN("Merging Tile Results");
				mergeTileResults();
			}

			for (auto& tile : m_Tiles) {
				tile.triangles.clear();
			}
			m_Triangles.clear();
		}
	}

	void TiledPipeline::initializeTiles()
	{
		ZoneScopedN("initializeTiles");
		int width = m_Framebuffer->getWidth();
		int height = m_Framebuffer->getHeight();
		int numTilesX = (width + TILE_SIZE - 1) / TILE_SIZE;
		int numTilesY = (height + TILE_SIZE - 1) / TILE_SIZE;

		m_Tiles.clear();
		m_Tiles.reserve(numTilesX * numTilesY);

		for (int ty = 0; ty < numTilesY; ++ty) {
			for (int tx = 0; tx < numTilesX; ++tx) {
				Tile tile;
				tile.startX = tx * TILE_SIZE;
				tile.startY = ty * TILE_SIZE;
				tile.endX = std::min((tx + 1) * TILE_SIZE, width);
				tile.endY = std::min((ty + 1) * TILE_SIZE, height);
				m_Tiles.push_back(tile);
			}
		}
	}

	bool TiledPipeline::triangleIntersectsTile(const Triangle& tri, const Tile& tile)
	{
		ZoneScopedN("triangleIntersectsTile");
		// Quick reject if bounding boxes don't overlap
		if (tri.maxX < tile.startX || tri.minX > tile.endX ||
			tri.maxY < tile.startY || tri.minY > tile.endY)
		{
			return false;
		}
		return true;
	}

	void TiledPipeline::binTrianglesToTiles()
	{
		ZoneScopedN("binTrianglesToTiles");
		int width = m_Framebuffer->getWidth();
		int height = m_Framebuffer->getHeight();

		int numTilesX = (width + TILE_SIZE - 1) / TILE_SIZE;

		for (auto& tri : m_Triangles) {
			//if (tri.isBackface()) continue;

		   // We figure out which tiles might overlap with this triangle
			int startTileX = std::max(0, (int)tri.minX / TILE_SIZE);
			int startTileY = std::max(0, (int)tri.minY / TILE_SIZE);
			int endTileX = std::min(numTilesX - 1, (int)tri.maxX / TILE_SIZE);
			int endTileY = std::min((int)m_Tiles.size() / numTilesX - 1,
				(int)tri.maxY / TILE_SIZE);

			for (int ty = startTileY; ty <= endTileY; ++ty) {
				for (int tx = startTileX; tx <= endTileX; ++tx) {
					int tileIdx = ty * numTilesX + tx;
					if (triangleIntersectsTile(tri, m_Tiles[tileIdx])) {
						m_Tiles[tileIdx].triangles.push_back(&tri);
					}
				}
			}
		}
	}

	void TiledPipeline::processTile(size_t tileIdx, ThreadLocalBuffers& buffers)
	{
		ZoneScopedN("processTile");
		Tile& tile = m_Tiles[tileIdx];
		if (tile.triangles.empty()) return;
		buffers.clear();

		VSTransformedTriangle vsTri;
		for (Triangle* tri : tile.triangles) {
			ZoneScopedN("Processing Triangle in Tile");
			m_Shader->material = tri->material;

			vsTri.vertices[0] = m_Shader->vertex(tri->vertices[0], 0);
			vsTri.vertices[1] = m_Shader->vertex(tri->vertices[1], 1);
			vsTri.vertices[2] = m_Shader->vertex(tri->vertices[2], 2);

			rasterizeTriangleInTile_EDGE(*tri, tile, buffers, vsTri);
		}

		// Store partial tile size in TileResult so we can do a correct merge
		TileResult& result = m_TileResults[tileIdx];
		result.startX = tile.startX;
		result.startY = tile.startY;
		result.endX = tile.endX;
		result.endY = tile.endY;

		// Copy the buffers out
		std::copy(buffers.colorBuffer.begin(), buffers.colorBuffer.end(),
			result.colorBuffer.begin());
		std::copy(buffers.depthBuffer.begin(), buffers.depthBuffer.end(),
			result.depthBuffer.begin());
	}

	void TiledPipeline::rasterizeTriangleInTile_AVX256(
		const Triangle& tri,
		const Tile& tile,
		ThreadLocalBuffers& buffers,
		const VSTransformedTriangle& vsTri)
	{
		ZoneScopedN("rasterizeTriangleInTile_AVX");

		// Compute bounding box clipped to the tile
		int startX = std::max(tile.startX, (int)std::floor(tri.minX));
		int startY = std::max(tile.startY, (int)std::floor(tri.minY));
		int endX = std::min(tile.endX, (int)std::ceil(tri.maxX));
		int endY = std::min(tile.endY, (int)std::ceil(tri.maxY));
		if (startX >= endX || startY >= endY) return;

		// Precompute area (signed)
		float area = (tri.screenPos[1].x - tri.screenPos[0].x) *
			(tri.screenPos[2].y - tri.screenPos[0].y) -
			(tri.screenPos[2].x - tri.screenPos[0].x) *
			(tri.screenPos[1].y - tri.screenPos[0].y);
		if (std::fabs(area) < 1e-12f) {
			return; // Degenerate triangle
		}
		float invArea = 1.0f / area;

		// Edge function coefficients (a, b, c)
		// Edge 0: from v1->v2
		float e0_a = tri.screenPos[1].y - tri.screenPos[2].y;
		float e0_b = tri.screenPos[2].x - tri.screenPos[1].x;
		float e0_c = tri.screenPos[1].x * tri.screenPos[2].y -
			tri.screenPos[2].x * tri.screenPos[1].y;

		// Edge 1: from v2->v0
		float e1_a = tri.screenPos[2].y - tri.screenPos[0].y;
		float e1_b = tri.screenPos[0].x - tri.screenPos[2].x;
		float e1_c = tri.screenPos[2].x * tri.screenPos[0].y -
			tri.screenPos[0].x * tri.screenPos[2].y;

		// Edge 2: from v0->v1
		float e2_a = tri.screenPos[0].y - tri.screenPos[1].y;
		float e2_b = tri.screenPos[1].x - tri.screenPos[0].x;
		float e2_c = tri.screenPos[0].x * tri.screenPos[1].y -
			tri.screenPos[1].x * tri.screenPos[0].y;

		// Prepare for AVX coverage
		__m256 vec_e0_a = _mm256_set1_ps(e0_a);
		__m256 vec_e1_a = _mm256_set1_ps(e1_a);
		__m256 vec_e2_a = _mm256_set1_ps(e2_a);
		__m256 vec_invArea = _mm256_set1_ps(invArea);

		// Pointers for depth/color
		float* depthBuffer = buffers.depthBuffer.data();
		uint8_t* colorBuffer = buffers.colorBuffer.data();

		// The offsets for 8 consecutive pixels (0..7) as floats
		__m256 pxOffsets = _mm256_setr_ps(0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f);

		for (int py = startY; py < endY; ++py)
		{
			float py_center = py + 0.5f;

			// Baseline coverage for the left edge at (startX + 0.5, py_center)
			float row_e0 = e0_a * (startX + 0.5f) + e0_b * py_center + e0_c;
			float row_e1 = e1_a * (startX + 0.5f) + e1_b * py_center + e1_c;
			float row_e2 = e2_a * (startX + 0.5f) + e2_b * py_center + e2_c;

			for (int px = startX; px < endX; px += 8)
			{
				// 1) Compute coverage in AVX registers for pixels [px .. px+7]
				__m256 cov_e0 = _mm256_add_ps(_mm256_set1_ps(row_e0),
					_mm256_mul_ps(vec_e0_a, pxOffsets));

				__m256 cov_e1 = _mm256_add_ps(_mm256_set1_ps(row_e1),
					_mm256_mul_ps(vec_e1_a, pxOffsets));

				__m256 cov_e2 = _mm256_add_ps(_mm256_set1_ps(row_e2),
					_mm256_mul_ps(vec_e2_a, pxOffsets));

				// 2) Compare coverage >= 0
				__m256 mask0 = _mm256_cmp_ps(cov_e0, _mm256_set1_ps(0.f), _CMP_GE_OQ);
				__m256 mask1 = _mm256_cmp_ps(cov_e1, _mm256_set1_ps(0.f), _CMP_GE_OQ);
				__m256 mask2 = _mm256_cmp_ps(cov_e2, _mm256_set1_ps(0.f), _CMP_GE_OQ);

				__m256 coverageMask = _mm256_and_ps(mask0, _mm256_and_ps(mask1, mask2));

				// 3) Convert to bitmask
				int maskBits = _mm256_movemask_ps(coverageMask);
				if (maskBits == 0) {
					// No pixels in this batch pass coverage
					// Move baseline coverage by 8 horizontally
					row_e0 += e0_a * 8.0f;
					row_e1 += e1_a * 8.0f;
					row_e2 += e2_a * 8.0f;
					continue;
				}

				// 4) Store coverage in arrays so we can access them per pixel
				float covE0[8], covE1[8], covE2[8];
				_mm256_storeu_ps(covE0, cov_e0);
				_mm256_storeu_ps(covE1, cov_e1);
				_mm256_storeu_ps(covE2, cov_e2);

				// 5) Process each pixel that passes coverage
				for (int i = 0; i < 8; ++i)
				{
					if (maskBits & (1 << i))
					{
						int curX = px + i;
						if (curX >= endX) break;  // boundary check

						float alpha = covE0[i] * invArea;
						float beta = covE1[i] * invArea;
						float gamma = covE2[i] * invArea;

						// Interpolate depth
						float z = tri.ndcZ[0] * alpha
							+ tri.ndcZ[1] * beta
							+ tri.ndcZ[2] * gamma;

						// Depth test
						int localX = curX - tile.startX;
						int localY = py - tile.startY;
						int idx = localY * TILE_SIZE + localX;

						if (z < depthBuffer[idx])
						{
							glm::vec4 colorOut;
							glm::vec3 bar(alpha, beta, gamma);
							bool discard = m_Shader->fragment(bar, colorOut, vsTri);

							if (!discard)
							{
								depthBuffer[idx] = z;
								int cIdx = idx * 4;
								colorBuffer[cIdx + 0] = (uint8_t)(std::clamp(colorOut.r, 0.0f, 1.0f) * 255.0f);
								colorBuffer[cIdx + 1] = (uint8_t)(std::clamp(colorOut.g, 0.0f, 1.0f) * 255.0f);
								colorBuffer[cIdx + 2] = (uint8_t)(std::clamp(colorOut.b, 0.0f, 1.0f) * 255.0f);
								colorBuffer[cIdx + 3] = (uint8_t)(std::clamp(colorOut.a, 0.0f, 1.0f) * 255.0f);
							}
						}
					}
				}

				// 6) Advance coverage row by 8 pixels
				row_e0 += e0_a * 8.0f;
				row_e1 += e1_a * 8.0f;
				row_e2 += e2_a * 8.0f;
			}
		}
	}

	void TiledPipeline::rasterizeTriangleInTile_AVX512(
		const Triangle& tri,
		const Tile& tile,
		ThreadLocalBuffers& buffers,
		const VSTransformedTriangle& vsTri)
	{
		ZoneScopedN("rasterizeTriangleInTile_AVX");

		// Compute bounding box clipped to the tile
		int startX = std::max(tile.startX, (int)std::floor(tri.minX));
		int startY = std::max(tile.startY, (int)std::floor(tri.minY));
		int endX = std::min(tile.endX, (int)std::ceil(tri.maxX));
		int endY = std::min(tile.endY, (int)std::ceil(tri.maxY));
		if (startX >= endX || startY >= endY) return;

		// Precompute area (signed)
		float area = (tri.screenPos[1].x - tri.screenPos[0].x) *
			(tri.screenPos[2].y - tri.screenPos[0].y) -
			(tri.screenPos[2].x - tri.screenPos[0].x) *
			(tri.screenPos[1].y - tri.screenPos[0].y);
		if (std::fabs(area) < 1e-12f) {
			return; // Degenerate triangle
		}
		float invArea = 1.0f / area;

		// Edge 0: from v1->v2
		float e0_a = tri.screenPos[1].y - tri.screenPos[2].y;
		float e0_b = tri.screenPos[2].x - tri.screenPos[1].x;
		float e0_c = tri.screenPos[1].x * tri.screenPos[2].y -
			tri.screenPos[2].x * tri.screenPos[1].y;

		// Edge 1: from v2->v0
		float e1_a = tri.screenPos[2].y - tri.screenPos[0].y;
		float e1_b = tri.screenPos[0].x - tri.screenPos[2].x;
		float e1_c = tri.screenPos[2].x * tri.screenPos[0].y -
			tri.screenPos[0].x * tri.screenPos[2].y;

		// Edge 2: from v0->v1
		float e2_a = tri.screenPos[0].y - tri.screenPos[1].y;
		float e2_b = tri.screenPos[1].x - tri.screenPos[0].x;
		float e2_c = tri.screenPos[0].x * tri.screenPos[1].y -
			tri.screenPos[1].x * tri.screenPos[0].y;

		// Prepare for AVX coverage
		__m512 vec_e0_a = _mm512_set1_ps(e0_a);
		__m512 vec_e1_a = _mm512_set1_ps(e1_a);
		__m512 vec_e2_a = _mm512_set1_ps(e2_a);
		__m512 vec_invArea = _mm512_set1_ps(invArea);

		// Pointers for depth/color
		float* depthBuffer = buffers.depthBuffer.data();
		uint8_t* colorBuffer = buffers.colorBuffer.data();

		// The offsets for 16 consecutive pixels (0..15) as floats
		__m512 pxOffsets = _mm512_setr_ps(0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f, 11.f, 12.f, 13.f, 14.f, 15.f);

		for (int py = startY; py < endY; ++py)
		{
			float py_center = py + 0.5f;

			// Baseline coverage for the left edge at (startX + 0.5, py_center)
			float row_e0 = e0_a * (startX + 0.5f) + e0_b * py_center + e0_c;
			float row_e1 = e1_a * (startX + 0.5f) + e1_b * py_center + e1_c;
			float row_e2 = e2_a * (startX + 0.5f) + e2_b * py_center + e2_c;

			for (int px = startX; px < endX; px += 8)
			{
				// 1) Compute coverage in AVX registers for pixels [px .. px+15]
				__m512 cov_e0 = _mm512_add_ps(_mm512_set1_ps(row_e0), _mm512_mul_ps(vec_e0_a, pxOffsets));

				__m512 cov_e1 = _mm512_add_ps(_mm512_set1_ps(row_e1),
					_mm512_mul_ps(vec_e1_a, pxOffsets));

				__m512 cov_e2 = _mm512_add_ps(_mm512_set1_ps(row_e2),
					_mm512_mul_ps(vec_e2_a, pxOffsets));

				// 2) Compare coverage >= 0
				__mmask16 mask0 = _mm512_cmp_ps_mask(cov_e0, _mm512_set1_ps(0.0f), _CMP_GE_OQ);
				__mmask16 mask1 = _mm512_cmp_ps_mask(cov_e1, _mm512_set1_ps(0.0f), _CMP_GE_OQ);
				__mmask16 mask2 = _mm512_cmp_ps_mask(cov_e2, _mm512_set1_ps(0.0f), _CMP_GE_OQ);

				__mmask16 coverageMask = mask0 & mask1 & mask2;

				// 3) Convert to bitmask
				int maskBits = _mm512_mask2int(coverageMask);
				if (maskBits == 0) {
					// No pixels in this batch pass coverage
					// Move baseline coverage by 8 horizontally
					row_e0 += e0_a * 16.0f;
					row_e1 += e1_a * 16.0f;
					row_e2 += e2_a * 16.0f;
					continue;
				}

				// 4) Store coverage in arrays so we can access them per pixel
				float covE0[16], covE1[16], covE2[16];
				_mm512_storeu_ps(covE0, cov_e0);
				_mm512_storeu_ps(covE1, cov_e1);
				_mm512_storeu_ps(covE2, cov_e2);

				// 5) Process each pixel that passes coverage
				for (int i = 0; i < 16; ++i)
				{
					if (maskBits & (1 << i))
					{
						int curX = px + i;
						if (curX >= endX) break;  // boundary check

						float alpha = covE0[i] * invArea;
						float beta = covE1[i] * invArea;
						float gamma = covE2[i] * invArea;

						glm::vec3 bar(alpha, beta, gamma);

						// Interpolate depth
						float z = tri.ndcZ[0] * alpha
							+ tri.ndcZ[1] * beta
							+ tri.ndcZ[2] * gamma;

						// Depth test
						int localX = curX - tile.startX;
						int localY = py - tile.startY;
						int idx = localY * TILE_SIZE + localX;

						if (z < depthBuffer[idx])
						{
							glm::vec4 colorOut;
							bool discard = m_Shader->fragment(bar,
								colorOut, vsTri);
							if (!discard)
							{
								depthBuffer[idx] = z;
								int cIdx = idx * 4;
								colorBuffer[cIdx + 0] = (uint8_t)(std::clamp(colorOut.r, 0.0f, 1.0f) * 255.0f);
								colorBuffer[cIdx + 1] = (uint8_t)(std::clamp(colorOut.g, 0.0f, 1.0f) * 255.0f);
								colorBuffer[cIdx + 2] = (uint8_t)(std::clamp(colorOut.b, 0.0f, 1.0f) * 255.0f);
								colorBuffer[cIdx + 3] = (uint8_t)(std::clamp(colorOut.a, 0.0f, 1.0f) * 255.0f);
							}
						}
					}
				}

				// 6) Advance coverage row by 8 pixels
				row_e0 += e0_a * 16.0f;
				row_e1 += e1_a * 16.0f;
				row_e2 += e2_a * 16.0f;
			}
		}
	}

	void TiledPipeline::rasterizeTriangleInTile_SSE(
		const Triangle& tri,
		const Tile& tile,
		ThreadLocalBuffers& buffers,
		const VSTransformedTriangle& vsTri)
	{
		ZoneScopedN("rasterizeTriangleInTile");

		// Compute bounding box clipped to the tile
		int startX = std::max(tile.startX, (int)std::floor(tri.minX));
		int startY = std::max(tile.startY, (int)std::floor(tri.minY));
		int endX = std::min(tile.endX, (int)std::ceil(tri.maxX));
		int endY = std::min(tile.endY, (int)std::ceil(tri.maxY));

		if (startX >= endX || startY >= endY) {
			return;
		}

		float x0 = tri.screenPos[0].x;
		float y0 = tri.screenPos[0].y;
		float x1 = tri.screenPos[1].x;
		float y1 = tri.screenPos[1].y;
		float x2 = tri.screenPos[2].x;
		float y2 = tri.screenPos[2].y;

		float area = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
		if (std::fabs(area) < 1e-12f) {
			return;  // Degenerate triangle
		}
		float invArea = 1.0f / area;

		// Edge 0: from vertex1->vertex2
		float e0_a = y1 - y2;
		float e0_b = x2 - x1;
		float e0_c = x1 * y2 - x2 * y1;

		// Edge 1: from vertex2->vertex0
		float e1_a = y2 - y0;
		float e1_b = x0 - x2;
		float e1_c = x2 * y0 - x0 * y2;

		// Edge 2: from vertex0->vertex1
		float e2_a = y0 - y1;
		float e2_b = x1 - x0;
		float e2_c = x0 * y1 - x1 * y0;

		// Local pointers for faster access
		float* depthBuffer = buffers.depthBuffer.data();
		uint8_t* colorBuffer = buffers.colorBuffer.data();

		__m128 vec_e0_a = _mm_set1_ps(e0_a);
		__m128 vec_e1_a = _mm_set1_ps(e1_a);
		__m128 vec_e2_a = _mm_set1_ps(e2_a);
		__m128 vec_invArea = _mm_set1_ps(invArea);

		for (int py = startY; py < endY; ++py)
		{
			// Center of the pixel row in the Y direction
			float py_center = py + 0.5f;

			// Edge function values at the left edge of the bounding box (the pixel center)
			float row_e0 = e0_a * (startX + 0.5f) + e0_b * py_center + e0_c;
			float row_e1 = e1_a * (startX + 0.5f) + e1_b * py_center + e1_c;
			float row_e2 = e2_a * (startX + 0.5f) + e2_b * py_center + e2_c;

			for (int px = startX; px < endX; px += 4)
			{
				// 1) Compute coverage for 4 consecutive pixels [px, px+1, px+2, px+3].
				//    We use (row_e0 + e0_a*i) to get each subpixel offset horizontally.
				__m128 pxOffsets = _mm_setr_ps(0.f, 1.f, 2.f, 3.f);

				// Coverage for edge 0
				__m128 coverage_e0 = _mm_add_ps(_mm_set1_ps(row_e0),
					_mm_mul_ps(vec_e0_a, pxOffsets));

				// Coverage for edge 1
				__m128 coverage_e1 = _mm_add_ps(_mm_set1_ps(row_e1),
					_mm_mul_ps(vec_e1_a, pxOffsets));

				// Coverage for edge 2
				__m128 coverage_e2 = _mm_add_ps(_mm_set1_ps(row_e2),
					_mm_mul_ps(vec_e2_a, pxOffsets));

				// 2) Compare coverage against 0 to see which pixels are inside
				__m128 mask0 = _mm_cmpge_ps(coverage_e0, _mm_set1_ps(0.f));
				__m128 mask1 = _mm_cmpge_ps(coverage_e1, _mm_set1_ps(0.f));
				__m128 mask2 = _mm_cmpge_ps(coverage_e2, _mm_set1_ps(0.f));
				__m128 coverageMask = _mm_and_ps(mask0, _mm_and_ps(mask1, mask2));

				// 3) Get bitmask for these 4 pixels
				int maskBits = _mm_movemask_ps(coverageMask);
				if (maskBits == 0)
				{
					// No pixels in this batch pass coverage
					// Move row coverage forward by 4 pixels
					row_e0 += e0_a * 4.0f;
					row_e1 += e1_a * 4.0f;
					row_e2 += e2_a * 4.0f;
					continue;
				}

				// 4) Store the coverage values in arrays, so we can access them per pixel
				float covE0[4], covE1[4], covE2[4];
				_mm_storeu_ps(covE0, coverage_e0);
				_mm_storeu_ps(covE1, coverage_e1);
				_mm_storeu_ps(covE2, coverage_e2);

				// 5) Process each pixel in the batch that is inside
				for (int i = 0; i < 4; ++i)
				{
					// If the i-th pixel is inside (bit i is set)
					if (maskBits & (1 << i))
					{
						int curX = px + i;
						if (curX >= endX) break;  // boundary check

						// Convert coverage to barycentric
						float alpha = covE0[i] * invArea;
						float beta = covE1[i] * invArea;
						float gamma = covE2[i] * invArea;

						glm::vec3 bar(alpha, beta, gamma);

						// Interpolate depth
						float z = tri.ndcZ[0] * alpha
							+ tri.ndcZ[1] * beta
							+ tri.ndcZ[2] * gamma;

						// Depth test
						int localX = curX - tile.startX;
						int localY = py - tile.startY;
						int idx = localY * TILE_SIZE + localX;

						if (z < depthBuffer[idx]) {
							glm::vec4 colorOut;
							bool discard = m_Shader->fragment(bar,
								colorOut, vsTri);
							if (!discard)
							{
								depthBuffer[idx] = z;
								int cIdx = idx * 4;
								colorBuffer[cIdx + 0] = (uint8_t)(std::clamp(colorOut.r, 0.0f, 1.0f) * 255.0f);
								colorBuffer[cIdx + 1] = (uint8_t)(std::clamp(colorOut.g, 0.0f, 1.0f) * 255.0f);
								colorBuffer[cIdx + 2] = (uint8_t)(std::clamp(colorOut.b, 0.0f, 1.0f) * 255.0f);
								colorBuffer[cIdx + 3] = (uint8_t)(std::clamp(colorOut.a, 0.0f, 1.0f) * 255.0f);
							}
						}
					}
				}

				row_e0 += e0_a * 4.0f;
				row_e1 += e1_a * 4.0f;
				row_e2 += e2_a * 4.0f;
			}
		}
	}

	void TiledPipeline::rasterizeTriangleInTile_EDGE(
		const Triangle& tri,
		const Tile& tile,
		ThreadLocalBuffers& buffers,
		const VSTransformedTriangle& vsTri)
	{
		ZoneScopedN("rasterizeTriangleInTile");

		// Compute bounding box clipped to the tile
		int startX = std::max(tile.startX, (int)std::floor(tri.minX));
		int startY = std::max(tile.startY, (int)std::floor(tri.minY));
		int endX = std::min(tile.endX, (int)std::ceil(tri.maxX));
		int endY = std::min(tile.endY, (int)std::ceil(tri.maxY));

		if (startX >= endX || startY >= endY) {
			return;
		}

		float x0 = tri.screenPos[0].x;
		float y0 = tri.screenPos[0].y;
		float x1 = tri.screenPos[1].x;
		float y1 = tri.screenPos[1].y;
		float x2 = tri.screenPos[2].x;
		float y2 = tri.screenPos[2].y;

		float area = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
		if (std::fabs(area) < 1e-12f) {
			return;  // Degenerate triangle
		}
		float invArea = 1.0f / area;

		// Edge 0: from vertex1->vertex2
		float e0_a = y1 - y2;
		float e0_b = x2 - x1;
		float e0_c = x1 * y2 - x2 * y1;

		// Edge 1: from vertex2->vertex0
		float e1_a = y2 - y0;
		float e1_b = x0 - x2;
		float e1_c = x2 * y0 - x0 * y2;

		// Edge 2: from vertex0->vertex1
		float e2_a = y0 - y1;
		float e2_b = x1 - x0;
		float e2_c = x0 * y1 - x1 * y0;

		// Local pointers for faster access
		float* depthBuffer = buffers.depthBuffer.data();
		uint8_t* colorBuffer = buffers.colorBuffer.data();

		for (int py = startY; py < endY; ++py)
		{
			float py_center = py + 0.5f;

			// Edge function values at (startX + 0.5, py_center)
			float row_e0 = e0_a * (startX + 0.5f) + e0_b * py_center + e0_c;
			float row_e1 = e1_a * (startX + 0.5f) + e1_b * py_center + e1_c;
			float row_e2 = e2_a * (startX + 0.5f) + e2_b * py_center + e2_c;

			for (int px = startX; px < endX; ++px)
			{
				float alpha = row_e0 * invArea;
				float beta = row_e1 * invArea;
				float gamma = row_e2 * invArea;

				if (alpha >= 0.f && beta >= 0.f && gamma >= 0.f)
				{
					glm::vec3 bar(alpha, beta, gamma);

					// Interpolate depth
					float z = tri.ndcZ[0] * alpha
						+ tri.ndcZ[1] * beta
						+ tri.ndcZ[2] * gamma;

					int localX = px - tile.startX;
					int localY = py - tile.startY;
					int idx = localY * TILE_SIZE + localX;
					if (z < depthBuffer[idx]) {
						glm::vec4 colorOut;
						bool discard = m_Shader->fragment(bar, colorOut, vsTri);
						if (!discard) {
							depthBuffer[idx] = z;
							int cIdx = idx * 4;
							colorBuffer[cIdx + 0] = (uint8_t)(std::clamp(colorOut.r, 0.0f, 1.0f) * 255.0f);
							colorBuffer[cIdx + 1] = (uint8_t)(std::clamp(colorOut.g, 0.0f, 1.0f) * 255.0f);
							colorBuffer[cIdx + 2] = (uint8_t)(std::clamp(colorOut.b, 0.0f, 1.0f) * 255.0f);
							colorBuffer[cIdx + 3] = (uint8_t)(std::clamp(colorOut.a, 0.0f, 1.0f) * 255.0f);
						}
					}
				}

				// Move to the next pixel horizontally
				row_e0 += e0_a;
				row_e1 += e1_a;
				row_e2 += e2_a;
			}
		}
	}

	void TiledPipeline::rasterizeTriangleInTile(const Triangle& tri,
		const Tile& tile,
		ThreadLocalBuffers& buffers, const VSTransformedTriangle& vsTri)
	{
		ZoneScopedN("rasterizeTriangleInTile");
		// Clamp bounding box to tile
		int startX = std::max(tile.startX, (int)std::floor(tri.minX));
		int startY = std::max(tile.startY, (int)std::floor(tri.minY));
		int endX = std::min(tile.endX, (int)std::ceil(tri.maxX));
		int endY = std::min(tile.endY, (int)std::ceil(tri.maxY));
		if (startX >= endX || startY >= endY) return;

		for (int py = startY; py < endY; ++py) {
			for (int px = startX; px < endX; ++px) {
				// Pixel center
				glm::vec2 p(px + 0.5f, py + 0.5f);

				glm::vec3 bcScreen = barycentric(
					glm::vec2(tri.screenPos[0].x, tri.screenPos[0].y),
					glm::vec2(tri.screenPos[1].x, tri.screenPos[1].y),
					glm::vec2(tri.screenPos[2].x, tri.screenPos[2].y),
					p
				);

				// If point is outside the triangle, skip
				if (bcScreen.x < 0.0f || bcScreen.y < 0.0f || bcScreen.z < 0.0f) {
					continue;
				}

				// Interpolate z/w and 1/w separately
				float z = tri.ndcZ[0] * bcScreen.x +
					tri.ndcZ[1] * bcScreen.y +
					tri.ndcZ[2] * bcScreen.z;

				// Convert to local tile coords
				int localX = px - tile.startX;
				int localY = py - tile.startY;
				int localIdx = localY * (TILE_SIZE)+localX;

				if (z < buffers.depthBuffer[localIdx])
				{
					glm::vec4 colorOut;

					// Fragment shader
					bool discard = m_Shader->fragment(bcScreen, colorOut, vsTri);
					if (!discard)
					{
						buffers.depthBuffer[localIdx] = z;
						int cIdx = localIdx * 4;
						buffers.colorBuffer[cIdx + 0] = static_cast<uint8_t>(
							std::clamp(colorOut.r, 0.0f, 1.0f) * 255.0f);
						buffers.colorBuffer[cIdx + 1] = static_cast<uint8_t>(
							std::clamp(colorOut.g, 0.0f, 1.0f) * 255.0f);
						buffers.colorBuffer[cIdx + 2] = static_cast<uint8_t>(
							std::clamp(colorOut.b, 0.0f, 1.0f) * 255.0f);
						buffers.colorBuffer[cIdx + 3] = static_cast<uint8_t>(
							std::clamp(colorOut.a, 0.0f, 1.0f) * 255.0f);
					}
				}
			}
		}
	}

	void TiledPipeline::mergeTileResults()
	{
		ZoneScopedN("mergeTileResults");
		for (size_t i = 0; i < m_TileResults.size(); ++i) {
			const auto& result = m_TileResults[i];
			int tileWidth = result.endX - result.startX;
			int tileHeight = result.endY - result.startY;

			// Merge each pixel of the tile result into the main framebuffer
			for (int y = 0; y < tileHeight; ++y)
			{
				int globalY = result.startY + y;
				for (int x = 0; x < tileWidth; ++x)
				{
					int globalX = result.startX + x;

					int localIdx = y * TILE_SIZE + x;
					float depth = result.depthBuffer[localIdx];
					if (depth < m_Framebuffer->getDepth(globalX, globalY))
					{
						int colorIdx = localIdx * 4;
						m_Framebuffer->setDepth(globalX, globalY, depth);
						m_Framebuffer->setPixel(globalX, globalY, {
							result.colorBuffer[colorIdx + 0],
							result.colorBuffer[colorIdx + 1],
							result.colorBuffer[colorIdx + 2],
							result.colorBuffer[colorIdx + 3]
							});
					}
				}
			}
		}
	}

	void TiledPipeline::workerThread()
	{
		ZoneScopedN("workerThread");
		ThreadLocalBuffers& buffers = t_buffers;

		while (true) {
			{
				std::unique_lock<std::mutex> lock(m_QueueMutex);
				m_WorkAvailable.wait(lock, [this]() {
					return m_HasWork || m_ShouldExit;
					});
				if (m_ShouldExit) return;
			}

			while (true) {
				int64_t tileIdx = m_NextTile.fetch_add(1, std::memory_order_relaxed);
				if (tileIdx >= m_TotalTiles) {
					break;
				}
				TracyPlot("Processing Tile", tileIdx);
				processTile(tileIdx, buffers);

				size_t doneCount = m_CompletedTiles.fetch_add(1) + 1;
				if (doneCount == m_TotalTiles) {
					m_WorkComplete.notify_one();
				}
			}
		}
	}
}