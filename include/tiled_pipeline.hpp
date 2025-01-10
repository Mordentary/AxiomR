#include "pipeline.hpp"
#include "thread_pool.hpp"
#include <mutex>
namespace AR
{
	class IShader;
	class Camera;
	class Framebuffer;
	struct VSTransformedTriangle;
	class Pipeline;

	struct Triangle {
		glm::vec2 screenPos[3];
		glm::vec3 ndcZ;
		std::array<Vertex, 3> vertices;
		float minX, minY, maxX, maxY;
		const Material* material;

		Triangle(const glm::vec4* clipSpace, const Vertex* verts, int width, int height, const Material* ptr);
		Triangle(const std::array<ClippedVertex, 3>& verts, int width, int height, const Material* mat);
		bool isBackface() const;
	};

	constexpr int TILE_SIZE = 16;
	constexpr size_t CACHE_LINE_SIZE = 64;
	constexpr int BATCH_SIZE = 8;

	template <size_t CurrentSize, size_t Alignment>
	struct Padding {
		static constexpr size_t padding_size = (CurrentSize % Alignment) == 0 ? 0 : Alignment - (CurrentSize % Alignment);
		std::array<char, padding_size> padding;

		Padding() = default;
	};
	template <size_t ColorBufferSize, size_t DepthBufferSize>
	struct InlinedBuffers {
		std::array<float, DepthBufferSize> depthBuffer;
		std::array<uint8_t, ColorBufferSize> colorBuffer; // RGBA for each pixel

		InlinedBuffers() {
			initialize();
		}

		void initialize() {
			std::memset(colorBuffer.data(), 0, colorBuffer.size());
			std::fill(depthBuffer.begin(), depthBuffer.end(), std::numeric_limits<float>::infinity());
		}

		void reset() {
			initialize();
		}
	};

	class TiledPipeline : public Pipeline {
	public:
		TiledPipeline(size_t threadsAvaibale, Camera* cam, Framebuffer* fb);
		~TiledPipeline();
		void drawMesh(const glm::mat4& modelMatrix, const Mesh& mesh) override;
	private:

		struct alignas(CACHE_LINE_SIZE) TileResult {
			int startX, startY, endX, endY;
			InlinedBuffers<TILE_SIZE* TILE_SIZE * 4, TILE_SIZE* TILE_SIZE> buffers;
			Padding<sizeof(InlinedBuffers<TILE_SIZE* TILE_SIZE * 4, TILE_SIZE* TILE_SIZE>) + 4 * sizeof(int), CACHE_LINE_SIZE> padding;

			TileResult()
				: startX(0), startY(0), endX(0), endY(0)
			{
			}

			void reset() {
				buffers.reset();
				startX = startY = endX = endY = 0;
			}
		};

		struct alignas(CACHE_LINE_SIZE) ThreadLocalBuffers {
			InlinedBuffers<TILE_SIZE* TILE_SIZE * 4, TILE_SIZE* TILE_SIZE> buffers;
			Padding<sizeof(InlinedBuffers<TILE_SIZE* TILE_SIZE * 4, TILE_SIZE* TILE_SIZE>), CACHE_LINE_SIZE> padding;

			ThreadLocalBuffers() {
				clear();
			}

			void clear() {
				buffers.reset();
			}
		};

		struct Tile {
			int startX, startY;
			int endX, endY;
			std::vector<Triangle*> triangles;
		};

		ThreadPool threadPool;

		std::vector<Tile> relevantTiles;
		std::atomic<size_t> currentBatchIndex;
		size_t totalBatches;

		// Helper methods
		void preprocessTiles(const std::vector<Tile>& allTiles);
		inline void createBatches() {
			totalBatches = (relevantTiles.size() + BATCH_SIZE - 1) / BATCH_SIZE;
		}
		// Tiles
		std::vector<Tile> m_Tiles;
		std::vector<TileResult> m_TileResults;

		// Triangle batch for current draw call
		std::vector<Triangle> m_Triangles;
		static thread_local ThreadLocalBuffers t_buffers;

		size_t m_NumThreads = 1;
		// Internal methods
		void initializeTiles();
		bool triangleIntersectsTile(const Triangle& tri, const Tile& tile);
		void binTrianglesToTiles();
		void processTile(size_t tileIdx, ThreadLocalBuffers& buffers);
		void rasterizeTriangleInTile(const Triangle& tri, const Tile& tile, ThreadLocalBuffers& buffers, const VSTransformedTriangle& vsTri);
		void rasterizeTriangleInTile_EDGE(const Triangle& tri, const Tile& tile, ThreadLocalBuffers& buffers, const VSTransformedTriangle& vsTri);
		void rasterizeTriangleInTile_AVX512(const Triangle& tri, const Tile& tile, ThreadLocalBuffers& buffers, const VSTransformedTriangle& vsTri);
		void rasterizeTriangleInTile_AVX256(const Triangle& tri, const Tile& tile, ThreadLocalBuffers& buffers, const VSTransformedTriangle& vsTri);
		void rasterizeTriangleInTile_SSE(const Triangle& tri, const Tile& tile, ThreadLocalBuffers& buffers, const VSTransformedTriangle& vsTri);

		void mergeTileResults();
	};
}