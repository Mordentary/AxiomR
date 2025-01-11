#include"tiled_pipeline.hpp"
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
#include <execution>
namespace AR
{
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
	thread_local TiledPipeline::ThreadLocalBuffers TiledPipeline::t_buffers{};

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

	Triangle::Triangle(const std::array<ClippedVertex, 3>& verts, int width, int height, const Material* mat)
		: material(mat)
	{
		for (size_t i = 0; i < verts.size(); ++i) {
			vertices[i] = verts[i].vertex;

			float w = clampW(verts[i].clipPos.w);
			float invW = 1.0f / w;

			ndcZ[i] = verts[i].clipPos.z * invW;

			screenPos[i].x = ((verts[i].clipPos.x * invW) + 1.0f) * 0.5f * width;
			screenPos[i].y = ((verts[i].clipPos.y * invW) + 1.0f) * 0.5f * height;
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

	TiledPipeline::TiledPipeline(size_t numThreads, Camera* cam, Framebuffer* fb)
		: m_ThreadPool(numThreads), m_CurrentBatchIndex(0), m_TotalBatches(0), m_NumThreads(numThreads), Pipeline(cam, fb)
	{
		initializeTiles();
	}
	TiledPipeline::~TiledPipeline()
	{
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

		const Vertex* vertices = mesh.getVertices().data();
		const std::vector<Face>& faces = mesh.getFaces();
		const auto& groups = mesh.getMaterialGroups();

		uint32_t framebufferWidth = m_Framebuffer->getWidth();
		uint32_t framebufferHeight = m_Framebuffer->getHeight();

		std::vector<Triangle> allClippedTriangles;
		allClippedTriangles.reserve(mesh.getFaces().size() * 2);

		{
			ZoneScopedN("ClippingAllTriangles");
			std::vector<std::future<void>> futures;

			for (const auto& group : groups)
			{
				ZoneScopedN("Clipping Material Group");
				const Material* material = mesh.getMaterial(group.materialName);

				// Distribute faces among threads
				int facesPerThread = (int)((group.faceCount + m_NumThreads - 1) / m_NumThreads);

				// For each group, create a per-thread clipped triangle array
				std::vector<std::vector<Triangle>> threadClipped(m_NumThreads);

				// Enqueue clipping jobs
				for (int t = 0; t < m_NumThreads; ++t)
				{
					futures.emplace_back(m_ThreadPool.enqueue([&, t]()
						{
							ZoneScopedN("ClipThread_Function");

							int start = (int)group.startIndex + t * facesPerThread;
							int end = std::min(start + facesPerThread, (int)(group.startIndex + group.faceCount));

							std::vector<Triangle> localTris;
							localTris.reserve(end - start);

							for (int i = start; i < end; ++i)
							{
								ZoneScopedN("ClipThread_Face");

								const auto& face = faces[i];
								assert(face.vertexIndices.size() == 3);

								const Vertex& v0 = vertices[face.vertexIndices[0]];
								const Vertex& v1 = vertices[face.vertexIndices[1]];
								const Vertex& v2 = vertices[face.vertexIndices[2]];

								{
									ZoneScopedN("ClipThread_TransformVertices");

									std::vector<ClippedVertex> clippedVertices;
									clippedVertices.reserve(8);
									clippedVertices.emplace_back(v0, mvp * glm::vec4(v0.position, 1.0f));
									clippedVertices.emplace_back(v1, mvp * glm::vec4(v1.position, 1.0f));
									clippedVertices.emplace_back(v2, mvp * glm::vec4(v2.position, 1.0f));

									{
										ZoneScopedN("ClipThread_ClipTriangle");
										clipTriangle(clippedVertices);
									}

									if (clippedVertices.size() >= 3)
									{
										ZoneScopedN("ClipThread_AssembleTriangles");
										for (size_t j = 1; j < clippedVertices.size() - 1; ++j)
										{
											ZoneScopedN("ClipThread_NewTriangle");
											const auto& cv0 = clippedVertices[0];
											const auto& cv1 = clippedVertices[j];
											const auto& cv2 = clippedVertices[j + 1];

											Triangle newTri{
												std::array<ClippedVertex, 3>{cv0, cv1, cv2},
												(int)framebufferWidth,
												(int)framebufferHeight,
												material
											};

											if (!newTri.isBackface())
											{
												localTris.push_back(std::move(newTri));
											}
										}
									}
								}

							}
								threadClipped[t] = std::move(localTris);
						}));
				}
				for (auto& f : futures) {
					f.get();
				}
				futures.clear();

				for (auto& tvec : threadClipped)
				{
					allClippedTriangles.insert(
						allClippedTriangles.end(),
						std::make_move_iterator(tvec.begin()),
						std::make_move_iterator(tvec.end())
					);
				}
			}
		}

		m_Triangles = std::move(allClippedTriangles);
		{
			ZoneScopedN("BinAllTriangles");
			binTrianglesToTiles();
		}

		preprocessTiles(m_Tiles);
		createBatches();
		m_TileResults.clear();
		m_TileResults.resize(m_RelevantTiles.size(), TileResult());
		m_CurrentBatchIndex.store(0, std::memory_order_relaxed);

		{
			ZoneScopedN("RasterizingTiles");

			std::vector<std::future<void>> futures;
			for (size_t t = 0; t < m_NumThreads; ++t) {
				futures.emplace_back(m_ThreadPool.enqueue([this]() {
					ZoneScopedN("TileWorker");
					while (true) {
						size_t batchIdx = m_CurrentBatchIndex.fetch_add(1, std::memory_order_relaxed);
						if (batchIdx >= m_TotalBatches)
							break;

						size_t startIdx = batchIdx * BATCH_SIZE;
						size_t endIdx = std::min(startIdx + BATCH_SIZE, m_RelevantTiles.size());

						for (size_t i = startIdx; i < endIdx; ++i) {
							processTile(i, t_buffers);
						}
					}
					}));
			}

			for (auto& f : futures) {
				f.get();
			}
		}

		{
			ZoneScopedN("MergingAllTiles");
			mergeTileResults();
		}

		m_RelevantTiles.clear();
		for (auto& tile : m_Tiles) {
			tile.triangles.clear();
		}
		m_Triangles.clear();
	}

	void TiledPipeline::preprocessTiles(const std::vector<Tile>& allTiles)
	{
		m_RelevantTiles.clear();
		m_RelevantTiles.reserve(allTiles.size());

		for (const auto& tile : allTiles) {
			if (!tile.triangles.empty()) {
				m_RelevantTiles.push_back(tile);
			}
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
		ZoneScoped;
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
		ZoneScoped;
		int width = m_Framebuffer->getWidth();
		int height = m_Framebuffer->getHeight();

		int numTilesX = (width + TILE_SIZE - 1) / TILE_SIZE;

		for (auto& tri : m_Triangles) {
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
		Tile& tile = m_RelevantTiles[tileIdx];
		if (tile.triangles.empty()) return;
		buffers.clear();

		VSTransformedTriangle vsTri;
		for (Triangle* tri : tile.triangles) {
			ZoneScopedN("Processing Triangle in Tile");
			m_Shader->material = tri->material;

			vsTri.vertices[0] = m_Shader->vertex(tri->vertices[0], 0);
			vsTri.vertices[1] = m_Shader->vertex(tri->vertices[1], 1);
			vsTri.vertices[2] = m_Shader->vertex(tri->vertices[2], 2);
			rasterizeTriangleInTile_AVX256(*tri, tile, buffers, vsTri);
		}

		// Store partial tile size in TileResult so we can do a correct merge
		TileResult& result = m_TileResults[tileIdx];
		result.startX = tile.startX;
		result.startY = tile.startY;
		result.endX = tile.endX;
		result.endY = tile.endY;

		// Copy the buffers out
		std::copy(buffers.buffers.colorBuffer.begin(), buffers.buffers.colorBuffer.end(),
			result.buffers.colorBuffer.begin());
		std::copy(buffers.buffers.depthBuffer.begin(), buffers.buffers.depthBuffer.end(),
			result.buffers.depthBuffer.begin());
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

		float x0 = tri.screenPos[0].x;
		float y0 = tri.screenPos[0].y;
		float x1 = tri.screenPos[1].x;
		float y1 = tri.screenPos[1].y;
		float x2 = tri.screenPos[2].x;
		float y2 = tri.screenPos[2].y;

		// Edge 0: from vertex1->vertex2
		float e0_c = x1 * y2 - x2 * y1;

		// Edge 1: from vertex2->vertex0
		float e1_c = x2 * y0 - x0 * y2;

		// Edge 2: from vertex0->vertex1
		float e2_c = x0 * y1 - x1 * y0;

		float area = e0_c + e1_c + e2_c;

		if (area >= 0 && area < 1.0E-12) return;

		// Edge 0: from vertex1->vertex2
		float e0_a = y1 - y2;
		float e0_b = x2 - x1;
		// Edge 1: from vertex2->vertex0
		float e1_a = y2 - y0;
		float e1_b = x0 - x2;
		// Edge 2: from vertex0->vertex1
		float e2_a = y0 - y1;
		float e2_b = x1 - x0;

		if ((area) < 0) {
			// Edge 1: inversion
			e0_a = -e0_a;
			e0_b = -e0_b;
			e0_c = -e0_c;
			// Edge 1: inversion
			e1_a = -e1_a;
			e1_b = -e1_b;
			e1_c = -e1_c;
			// Edge 2: inversion
			e2_a = -e2_a;
			e2_b = -e2_b;
			e2_c = -e2_c;
			area = -area;
		}
		float invArea = 1.0f / area;

		// Prepare for AVX coverage
		__m256 vec_e0_a = _mm256_set1_ps(e0_a);
		__m256 vec_e1_a = _mm256_set1_ps(e1_a);
		__m256 vec_e2_a = _mm256_set1_ps(e2_a);
		__m256 vec_invArea = _mm256_set1_ps(invArea);

		// Pointers for depth/color
		float* depthBuffer = buffers.buffers.depthBuffer.data();
		uint8_t* colorBuffer = buffers.buffers.colorBuffer.data();

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

		float x0 = tri.screenPos[0].x;
		float y0 = tri.screenPos[0].y;
		float x1 = tri.screenPos[1].x;
		float y1 = tri.screenPos[1].y;
		float x2 = tri.screenPos[2].x;
		float y2 = tri.screenPos[2].y;

		// Edge 0: from vertex1->vertex2
		float e0_c = x1 * y2 - x2 * y1;

		// Edge 1: from vertex2->vertex0
		float e1_c = x2 * y0 - x0 * y2;

		// Edge 2: from vertex0->vertex1
		float e2_c = x0 * y1 - x1 * y0;

		float area = e0_c + e1_c + e2_c;

		if (area >= 0 && area < 1.0E-12) return;

		// Edge 0: from vertex1->vertex2
		float e0_a = y1 - y2;
		float e0_b = x2 - x1;
		// Edge 1: from vertex2->vertex0
		float e1_a = y2 - y0;
		float e1_b = x0 - x2;
		// Edge 2: from vertex0->vertex1
		float e2_a = y0 - y1;
		float e2_b = x1 - x0;

		if ((area) < 0) {
			// Edge 1: inversion
			e0_a = -e0_a;
			e0_b = -e0_b;
			e0_c = -e0_c;
			// Edge 1: inversion
			e1_a = -e1_a;
			e1_b = -e1_b;
			e1_c = -e1_c;
			// Edge 2: inversion
			e2_a = -e2_a;
			e2_b = -e2_b;
			e2_c = -e2_c;
			area = -area;
		}
		float invArea = 1.0f / area;
		// Prepare for AVX coverage
		__m512 vec_e0_a = _mm512_set1_ps(e0_a);
		__m512 vec_e1_a = _mm512_set1_ps(e1_a);
		__m512 vec_e2_a = _mm512_set1_ps(e2_a);
		__m512 vec_invArea = _mm512_set1_ps(invArea);

		// Pointers for depth/color
		float* depthBuffer = buffers.buffers.depthBuffer.data();
		uint8_t* colorBuffer = buffers.buffers.colorBuffer.data();

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

		// Edge 0: from vertex1->vertex2
		float e0_c = x1 * y2 - x2 * y1;

		// Edge 1: from vertex2->vertex0
		float e1_c = x2 * y0 - x0 * y2;

		// Edge 2: from vertex0->vertex1
		float e2_c = x0 * y1 - x1 * y0;

		float area = e0_c + e1_c + e2_c;

		if (area >= 0 && area < 1.0E-12) return;

		// Edge 0: from vertex1->vertex2
		float e0_a = y1 - y2;
		float e0_b = x2 - x1;
		// Edge 1: from vertex2->vertex0
		float e1_a = y2 - y0;
		float e1_b = x0 - x2;
		// Edge 2: from vertex0->vertex1
		float e2_a = y0 - y1;
		float e2_b = x1 - x0;

		if ((area) < 0) {
			// Edge 1: inversion
			e0_a = -e0_a;
			e0_b = -e0_b;
			e0_c = -e0_c;
			// Edge 1: inversion
			e1_a = -e1_a;
			e1_b = -e1_b;
			e1_c = -e1_c;
			// Edge 2: inversion
			e2_a = -e2_a;
			e2_b = -e2_b;
			e2_c = -e2_c;
			area = -area;
		}
		float invArea = 1.0f / area;

		// Local pointers for faster access
		float* depthBuffer = buffers.buffers.depthBuffer.data();
		uint8_t* colorBuffer = buffers.buffers.colorBuffer.data();

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

		// Edge 0: from vertex1->vertex2
		float e0_c = x1 * y2 - x2 * y1;

		// Edge 1: from vertex2->vertex0
		float e1_c = x2 * y0 - x0 * y2;

		// Edge 2: from vertex0->vertex1
		float e2_c = x0 * y1 - x1 * y0;

		float area = e0_c + e1_c + e2_c;

		if (area >= 0 && area < 1.0E-12) return;

		// Edge 0: from vertex1->vertex2
		float e0_a = y1 - y2;
		float e0_b = x2 - x1;
		// Edge 1: from vertex2->vertex0
		float e1_a = y2 - y0;
		float e1_b = x0 - x2;
		// Edge 2: from vertex0->vertex1
		float e2_a = y0 - y1;
		float e2_b = x1 - x0;

		if ((area) < 0) {
			// Edge 1: inversion
			e0_a = -e0_a;
			e0_b = -e0_b;
			e0_c = -e0_c;
			// Edge 1: inversion
			e1_a = -e1_a;
			e1_b = -e1_b;
			e1_c = -e1_c;
			// Edge 2: inversion
			e2_a = -e2_a;
			e2_b = -e2_b;
			e2_c = -e2_c;
			area = -area;
		}
		float invArea = 1.0f / area;

		// Local pointers for faster access
		float* depthBuffer = buffers.buffers.depthBuffer.data();
		uint8_t* colorBuffer = buffers.buffers.colorBuffer.data();

		for (int py = startY; py < endY; ++py)
		{
			float py_center = py + 0.5f;

			// Edge function values at (startX + 0.5, py_center)
			float row_e0 = e0_a * (startX + 0.5f) + e0_b * py_center + e0_c;
			float row_e1 = e1_a * (startX + 0.5f) + e1_b * py_center + e1_c;
			float row_e2 = e2_a * (startX + 0.5f) + e2_b * py_center + e2_c;

			//float row_e0 = e0_a * (startX + 0.5f) + e0_b * py_center + e0_c;
			//float row_e1 = e1_a * (startX + 0.5f) + e1_b * py_center + e1_c;
			//float row_e2 = e2_a * (startX + 0.5f) + e2_b * py_center + e2_c;

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

				if (z < buffers.buffers.depthBuffer[localIdx])
				{
					glm::vec4 colorOut;

					// Fragment shader
					bool discard = m_Shader->fragment(bcScreen, colorOut, vsTri);
					if (!discard)
					{
						buffers.buffers.depthBuffer[localIdx] = z;
						int cIdx = localIdx * 4;
						buffers.buffers.colorBuffer[cIdx + 0] = static_cast<uint8_t>(
							std::clamp(colorOut.r, 0.0f, 1.0f) * 255.0f);
						buffers.buffers.colorBuffer[cIdx + 1] = static_cast<uint8_t>(
							std::clamp(colorOut.g, 0.0f, 1.0f) * 255.0f);
						buffers.buffers.colorBuffer[cIdx + 2] = static_cast<uint8_t>(
							std::clamp(colorOut.b, 0.0f, 1.0f) * 255.0f);
						buffers.buffers.colorBuffer[cIdx + 3] = static_cast<uint8_t>(
							std::clamp(colorOut.a, 0.0f, 1.0f) * 255.0f);
					}
				}
			}
		}
	}

	void TiledPipeline::mergeTileResults()
	{
		ZoneScopedN("mergeTileResults");

		float* mainDepth = m_Framebuffer->getDepthData();
		uint8_t* mainColor = m_Framebuffer->getColorData();

		int framebufferWidth = m_Framebuffer->getWidth();

		for (size_t i = 0; i < m_TileResults.size(); ++i) {
			const auto& result = m_TileResults[i];
			int tileWidth = result.endX - result.startX;
			int tileHeight = result.endY - result.startY;

			const float* depthData = result.buffers.depthBuffer.data();
			const uint8_t* colorData = result.buffers.colorBuffer.data();

			for (int y = 0; y < tileHeight; ++y)
			{
				int globalY = result.startY + y;
				for (int x = 0; x < tileWidth; x += 8) // Process 8 pixels at a time
				{
					int remaining = tileWidth - x >= 8 ? 8 : tileWidth - x;
					__m256 tileDepth = _mm256_loadu_ps(&depthData[y * TILE_SIZE + x]);
					__m256 mainDepthVals = _mm256_loadu_ps(&mainDepth[globalY * framebufferWidth + result.startX + x]);
					__m256 cmp = _mm256_cmp_ps(tileDepth, mainDepthVals, _CMP_LT_OS);
					int mask = _mm256_movemask_ps(cmp);

					if (mask != 0)
					{
						__m256 updatedDepth = _mm256_blendv_ps(mainDepthVals, tileDepth, cmp);
						_mm256_storeu_ps(&mainDepth[globalY * framebufferWidth + result.startX + x], updatedDepth);

						for (int bit = 0; bit < 8; ++bit)
						{
							if (mask & (1 << bit))
							{
								int pixelX = x + bit;
								if (pixelX >= tileWidth) break;

								int localIdx = y * TILE_SIZE + pixelX;
								int globalIdx = globalY * framebufferWidth + result.startX + pixelX;

								int tileColorIdx = localIdx * 4;
								int finalColorIdx = globalIdx * 4;

								mainColor[finalColorIdx + 0] = colorData[tileColorIdx + 2];
								mainColor[finalColorIdx + 1] = colorData[tileColorIdx + 1];
								mainColor[finalColorIdx + 2] = colorData[tileColorIdx + 0];
								mainColor[finalColorIdx + 3] = colorData[tileColorIdx + 3];
							}
						}
					}
				}
			}
		}
	}
}