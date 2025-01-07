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
		void clipTriangle(std::vector<std::pair<Vertex, glm::vec4>>& clippedVertices);
		void clipAgainstPlane(std::vector<std::pair<Vertex, glm::vec4>>& poly, int plane);
		inline bool insidePlane(const glm::vec4& v, int plane);
		inline float intersectPlane(const glm::vec4& v1, const glm::vec4& v2, int plane);
		inline std::pair<Vertex, glm::vec4> interpolateVertices(std::pair<Vertex, glm::vec4> v0, std::pair<Vertex, glm::vec4> v1, float t_Point);
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
		VSTransformedTriangle vsOutTriangle;

		Triangle(const glm::vec4* clipSpace, const Vertex* verts, int width, int height, const Material* ptr);
		Triangle(const std::array<std::pair<Vertex, glm::vec4>, 3>& verts, int width, int height, const Material* mat);
		bool isBackface() const;
	};

	class TiledPipeline : public Pipeline {
	private:
		static constexpr int TILE_SIZE = 32;

		struct TileResult {
			std::vector<uint8_t> colorBuffer;
			std::vector<float> depthBuffer;
			// Actual region we cover in final merge:
			int startX, startY, endX, endY;

			TileResult(int tileSize)
				: colorBuffer(tileSize* tileSize * 4, 0),
				depthBuffer(tileSize* tileSize, std::numeric_limits<float>::infinity()),
				startX(0), startY(0), endX(0), endY(0)
			{
			}
		};

		struct Tile {
			int startX, startY;
			int endX, endY;
			std::vector<Triangle*> triangles;
		};

		struct ThreadLocalBuffers {
			alignas(64) std::vector<uint8_t> colorBuffer;
			alignas(64) std::vector<float> depthBuffer;

			ThreadLocalBuffers(int tileSize)
				: colorBuffer(tileSize* tileSize * 4, 0),
				depthBuffer(tileSize* tileSize, std::numeric_limits<float>::infinity())
			{
			}

			void clear() {
				std::fill(colorBuffer.begin(), colorBuffer.end(), 0);
				std::fill(depthBuffer.begin(), depthBuffer.end(), std::numeric_limits<float>::infinity());
			}
		};

		// Thread pool
		std::vector<std::thread> m_Workers;
		// Tiles
		std::vector<Tile> m_Tiles;
		std::vector<TileResult> m_TileResults;

		// Triangle batch for current draw call
		std::vector<Triangle> m_Triangles;

		// Thread local storage
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
		void rasterizeTriangleInTile(const Triangle& tri, const Tile& tile, ThreadLocalBuffers& buffers);
		static float edgeFunction(const glm::vec2& a, const glm::vec2& b, const glm::vec2& c);
		void mergeTileResults();
		void mergeTileToFramebuffer(const Tile& tile, const ThreadLocalBuffers& buffers);
		void workerThread();

	public:
		TiledPipeline(int numThreads = std::thread::hardware_concurrency());
		~TiledPipeline();
		void drawMesh(const glm::mat4& modelMatrix, const Mesh& mesh) override;

		// Additional helper structures
		//struct EdgeEquation {
		//	float a, b, c;
		//	EdgeEquation(const glm::vec2& v0, const glm::vec2& v1);
		//	float evaluate(float x, float y) const;
		//};

		//struct TriangleSetup {
		//	EdgeEquation edges[3];
		//	float area;
		//	float oneOverArea;

		//	TriangleSetup(const Triangle& tri);
		//	bool isBackface() const;
		//	glm::vec3 computeBarycentrics(float x, float y) const;
		//};

		void optimizedRasterizeTriangleInTile(const Triangle& tri, const Tile& tile, ThreadLocalBuffers& buffers);
	};
}