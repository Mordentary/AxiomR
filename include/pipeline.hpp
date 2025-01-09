#pragma once
#include <algorithm>
#include <limits>
#include "IShader.hpp"
#include <mutex>
namespace AR
{
	class IShader;
	class Camera;
	class Framebuffer;
	struct VSTransformedTriangle;
	struct ClippedVertex {
		Vertex vertex;
		glm::vec4 clipPos;
	};
	class Pipeline {
		friend IShader;
	public:
		Pipeline();
		void setShader(IShader* shader);
		void setCamera(Camera* cam);
		void setFramebuffer(Framebuffer* fb);
		virtual void drawMesh(const glm::mat4& modelMatrix, const Mesh& mesh);
		glm::mat4 getViewportMat();
	protected:
		IShader* m_Shader = nullptr;
		Framebuffer* m_Framebuffer = nullptr;
		const Camera* m_Camera = nullptr;
		std::mutex m_Mutex;
		glm::vec3 barycentric(const glm::vec2& A, const glm::vec2& B, const glm::vec2& C, const glm::vec2& P) const;
		void rasterizeTriangle(const glm::vec4 clip[3]);
		void clipTriangle(std::vector<ClippedVertex>& clippedVertices);
		void clipAgainstPlane(std::vector<ClippedVertex>& poly, int plane, std::vector<ClippedVertex>& tempOut);
		inline bool insidePlane(const glm::vec4& v, int plane);
		inline float intersectPlane(const glm::vec4& v1, const glm::vec4& v2, int plane);
		inline ClippedVertex interpolateVertices(const ClippedVertex& v0, const ClippedVertex& v1, float t_Point);
		//void rasterizeTriangleInTile(const glm::vec4 clip[3], size_t threadId,
			//std::vector<std::vector<std::vector<float>>>& depthBuffers,
			//std::vector<std::vector<std::vector<glm::vec4>>>& colorBuffers,
			//int tileX, int tileY, int tileIndex);
	};

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

	class TiledPipeline : public Pipeline {
	private:
		static constexpr int TILE_SIZE = 16;
		static constexpr size_t CACHE_LINE_SIZE = 64;
		template <size_t CurrentSize, size_t Alignment>
		struct Padding {
			static constexpr size_t padding_size = (CurrentSize % Alignment) == 0 ? 0 : Alignment - (CurrentSize % Alignment);
			std::array<char, padding_size> padding;

			Padding() = default;
		};
		template <size_t ColorBufferSize, size_t DepthBufferSize>
		struct InlinedBuffers {
			std::array<uint8_t, ColorBufferSize> colorBuffer; // RGBA for each pixel
			std::array<float, DepthBufferSize> depthBuffer;

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

		std::vector<std::thread> m_Workers;
		// Tiles
		std::vector<Tile> m_Tiles;
		std::vector<TileResult> m_TileResults;

		// Triangle batch for current draw call
		std::vector<Triangle> m_Triangles;
		static thread_local ThreadLocalBuffers t_buffers;

		// Sync
		std::mutex m_QueueMutex;
		std::condition_variable m_WorkAvailable;
		std::condition_variable m_WorkComplete;
		std::atomic<bool> m_ShouldExit{ false };
		std::atomic<bool> m_HasWork{ false };
		std::atomic<size_t> m_CompletedTiles{ 0 };
		std::atomic<size_t> m_TotalTiles{ 0 };
		std::atomic<size_t> m_NextTile{ 0 };

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

		static float edgeFunction(const glm::vec2& a, const glm::vec2& b, const glm::vec2& c);
		void mergeTileResults();
		void mergeTileToFramebuffer(const Tile& tile, const ThreadLocalBuffers& buffers);
		void workerThread();

	public:
		TiledPipeline(int numThreads = std::thread::hardware_concurrency());
		~TiledPipeline();
		void drawMesh(const glm::mat4& modelMatrix, const Mesh& mesh) override;
	};
}