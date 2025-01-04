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
#include <vector>
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
			//	for (auto& tri : clippedTriangles) {
			//		glm::vec4 clipCoords[3];
			//		clipCoords[0] = m_Shader->vertex(tri[0].first, 0);
			//		clipCoords[1] = m_Shader->vertex(tri[1].first, 1);
			//		clipCoords[2] = m_Shader->vertex(tri[2].first, 2);

			//		rasterizeTriangle(clipCoords);
			//	}
			//}
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

	// Thread-local buffers definition
	thread_local TiledPipeline::ThreadLocalBuffers TiledPipeline::t_buffers(TILE_SIZE);

	// ---------------- Triangle methods ----------------
	Triangle::Triangle(const glm::vec4* clipSpace, const Vertex* verts,
		int width, int height, const Material* mat)
		: material(mat)
	{
		for (int i = 0; i < 3; i++) {
			vertices[i] = verts[i];
			clipPos[i] = clipSpace[i];

			float w = clipSpace[i].w;
			if (std::abs(w) < 1e-6f) {
				w = std::copysign(1e-6f, w);
			}
			if (std::abs(w) < 1e-6f) w = (w < 0.0f) ? -1e-6f : 1e-6f;

			float zNDC = clipSpace[i].z / w;
			ndcZ[i] = zNDC;

			screenPos[i].x = ((clipSpace[i].x / w) + 1.0f) * 0.5f * width;
			screenPos[i].y = ((clipSpace[i].y / w) + 1.0f) * 0.5f * height;
		}
		minX = std::min({ screenPos[0].x, screenPos[1].x, screenPos[2].x });
		minY = std::min({ screenPos[0].y, screenPos[1].y, screenPos[2].y });
		maxX = std::max({ screenPos[0].x, screenPos[1].x, screenPos[2].x });
		maxY = std::max({ screenPos[0].y, screenPos[1].y, screenPos[2].y });
	}

	bool Triangle::isBackface() const {
		// If you do 2D cross culling, do it here:
		// float area = ...
		// return (area <= 0.0f);
		return false;
	}

	// ---------------- TiledPipeline methods ----------------

	TiledPipeline::TiledPipeline(int numThreads)
	{
		// Spawn threads
		for (int i = 0; i < numThreads; ++i) {
			m_Workers.emplace_back([this]() { workerThread(); });
		}
	}

	TiledPipeline::~TiledPipeline()
	{
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
		if (!m_Shader || !m_Camera || !m_Framebuffer) return;

		// Set up transforms
		const glm::mat4& viewProj = m_Camera->getViewProjectionMatrix();
		glm::mat4 mvp = viewProj * modelMatrix;

		m_Shader->model = modelMatrix;
		m_Shader->viewProj = viewProj;
		m_Shader->mvp = mvp;
		m_Shader->viewportMat = m_Camera->getViewportMatrix();
		m_Shader->cameraPosition = m_Camera->getPosition();

		// Prepare tiles
		initializeTiles();
		m_TileResults.clear();
		m_TileResults.resize(m_Tiles.size(), TileResult(TILE_SIZE));

		// Build triangle list
		m_Triangles.clear();
		const auto& vertices = mesh.getVertices();
		const auto& faces = mesh.getFaces();
		const auto& groups = mesh.getMaterialGroups();

		for (auto& group : groups)
		{
			const Material* groupMat = mesh.getMaterial(group.materialName);
			for (int faceIdx = group.startIndex;
				faceIdx < group.startIndex + group.faceCount; ++faceIdx)
			{
				glm::vec4 clipSpace[3];
				Vertex triVerts[3];
				const Face& face = faces[faceIdx];
				for (int i = 0; i < 3; i++) {
					triVerts[i] = vertices[face.vertexIndices[i]];
					clipSpace[i] = mvp * glm::vec4(triVerts[i].position, 1.0f);
				}

				m_Triangles.emplace_back(clipSpace, triVerts,
					m_Framebuffer->getWidth(),
					m_Framebuffer->getHeight(),
					groupMat);
			}
		}

		// Bin to tiles
		binTrianglesToTiles();

		// Start worker threads
		m_NextTile = 0;
		m_CompletedTiles = 0;
		m_TotalTiles = m_Tiles.size();

		{
			std::unique_lock<std::mutex> lock(m_QueueMutex);
			m_HasWork = true;
		}
		m_WorkAvailable.notify_all();

		// Wait for finish
		{
			std::unique_lock<std::mutex> lock(m_QueueMutex);
			m_WorkComplete.wait(lock, [this]() {
				return m_CompletedTiles == m_TotalTiles;
				});
		}

		// Reset
		m_HasWork = false;

		// Merge results into the main framebuffer
		mergeTileResults();
	}

	void TiledPipeline::initializeTiles()
	{
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
		int width = m_Framebuffer->getWidth();
		int height = m_Framebuffer->getHeight();

		int numTilesX = (width + TILE_SIZE - 1) / TILE_SIZE;

		for (auto& tri : m_Triangles) {
			// if (tri.isBackface()) continue; // optional backface cull

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
		buffers.clear();
		Tile& tile = m_Tiles[tileIdx];

		for (Triangle* tri : tile.triangles) {
			m_Shader->material = tri->material;

			// Vertex shader outputs
			VertexOutput vOut0 = m_Shader->vertex(tri->vertices[0], 0);
			VertexOutput vOut1 = m_Shader->vertex(tri->vertices[1], 1);
			VertexOutput vOut2 = m_Shader->vertex(tri->vertices[2], 2);
			VSTransformedTriangle vsTri{ vOut0, vOut1, vOut2 };
			tri->vsOutTriangle = std::move(vsTri);

			rasterizeTriangleInTile(*tri, tile, buffers);
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

	// A simple edge function
	float TiledPipeline::edgeFunction(const glm::vec2& a,
		const glm::vec2& b,
		const glm::vec2& c)
	{
		return (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x);
	}

	void TiledPipeline::rasterizeTriangleInTile(const Triangle& tri,
		const Tile& tile,
		ThreadLocalBuffers& buffers)
	{
		// Clamp bounding box to tile
		int startX = std::max(tile.startX, (int)std::floor(tri.minX));
		int startY = std::max(tile.startY, (int)std::floor(tri.minY));
		int endX = std::min(tile.endX, (int)std::ceil(tri.maxX));
		int endY = std::min(tile.endY, (int)std::ceil(tri.maxY));
		if (startX >= endX || startY >= endY) return;

		float area = edgeFunction(tri.screenPos[0], tri.screenPos[1], tri.screenPos[2]);
		if (std::fabs(area) < 1e-6f) return;
		float invArea = 1.0f / area;

		for (int py = startY; py < endY; ++py) {
			for (int px = startX; px < endX; ++px) {
				// Pixel center
				glm::vec2 p(px + 0.5f, py + 0.5f);

				float w0 = edgeFunction(tri.screenPos[1], tri.screenPos[2], p);
				float w1 = edgeFunction(tri.screenPos[2], tri.screenPos[0], p);
				float w2 = edgeFunction(tri.screenPos[0], tri.screenPos[1], p);

				if (w0 >= 0.f && w1 >= 0.f && w2 >= 0.f) {
					w0 *= invArea;
					w1 *= invArea;
					w2 *= invArea;

					// Interpolate depth
					float z = tri.ndcZ[0] * w0 +
						tri.ndcZ[1] * w1 +
						tri.ndcZ[2] * w2;

					// Convert to local tile coords
					int localX = px - tile.startX;
					int localY = py - tile.startY;
					int localIdx = localY * TILE_SIZE + localX;

					if (z < buffers.depthBuffer[localIdx]) {
						glm::vec4 colorOut;
						glm::vec3 bary(w0, w1, w2);

						// Fragment shader
						bool discard = m_Shader->fragment(bary, colorOut, tri.vsOutTriangle);
						if (!discard) {
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
	}

	void TiledPipeline::mergeTileResults()
	{
		for (size_t i = 0; i < m_TileResults.size(); ++i) {
			const auto& result = m_TileResults[i];
			int tileWidth = result.endX - result.startX;
			int tileHeight = result.endY - result.startY;

			// Merge each pixel of the tile result into the main framebuffer
			for (int y = 0; y < tileHeight; ++y) {
				int globalY = result.startY + y;
				for (int x = 0; x < tileWidth; ++x) {
					int globalX = result.startX + x;

					int localIdx = y * TILE_SIZE + x; // we used a fixed size tile buffer
					float depth = result.depthBuffer[localIdx];
					if (depth <= m_Framebuffer->getDepth(globalX, globalY)) {
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
				size_t tileIdx = m_NextTile.fetch_add(1, std::memory_order_relaxed);
				if (tileIdx >= m_TotalTiles) {
					break;
				}
				buffers.clear();
				processTile(tileIdx, buffers);

				// Mark tile done
				size_t doneCount = m_CompletedTiles.fetch_add(1) + 1;
				if (doneCount == m_TotalTiles) {
					m_WorkComplete.notify_one();
				}
			}
		}
	}
}