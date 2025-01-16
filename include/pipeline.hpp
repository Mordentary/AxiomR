#pragma once
#include <algorithm>
#include<array>
#include <limits>
#include "IShader.hpp"
namespace AR
{
	class IShader;
	class Camera;
	class Framebuffer;
	struct VSTransformedTriangle;
	struct Triangle;
	struct ClippedVertex {
		Vertex vertex;
		glm::vec4 clipPos;
	};
	constexpr int MAX_CLIPPED_VERTS = 24;
	class Pipeline {
		friend IShader;
	public:
		Pipeline(Camera* cam, Framebuffer* fb);
		void setShader(IShader* shader);
		void setCamera(Camera* cam);
		void setFramebuffer(Framebuffer* fb);
		virtual void drawMesh(const glm::mat4& modelMatrix, const Mesh& mesh);
		glm::mat4 getViewportMat();
	protected:
		IShader* m_Shader = nullptr;
		Framebuffer* m_Framebuffer = nullptr;
		const Camera* m_Camera = nullptr;
		glm::vec3 barycentric(const glm::vec2& A, const glm::vec2& B, const glm::vec2& C, const glm::vec2& P) const;
		void rasterizeTriangle(const glm::vec4 clip[3]);

		static int clipTriangleSinglePlane(
			int plane,
			ClippedVertex& v0,
			ClippedVertex& v1,
			ClippedVertex& v2,
			ClippedVertex& v3);

		static void clipTriangle(size_t& arrayIndex,
			std::array<ClippedVertex, MAX_CLIPPED_VERTS>& outTris);

		static void clipAgainstPlane(std::pmr::vector<ClippedVertex>& poly, int plane, std::pmr::vector<ClippedVertex>& tempOut);
		static inline bool insidePlane(const glm::vec4& v, int plane);
		static inline float intersectPlane(const glm::vec4& v1, const glm::vec4& v2, int plane);
		static inline ClippedVertex interpolateVertices(const ClippedVertex& v0, const ClippedVertex& v1, float t_Point);
	};
}