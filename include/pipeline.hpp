#pragma once
#include <algorithm>
#include <limits>
#include "IShader.hpp"
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
		void clipTriangle(std::vector<ClippedVertex>& clippedVertices);
		void clipAgainstPlane(std::vector<ClippedVertex>& poly, int plane, std::vector<ClippedVertex>& tempOut);
		inline bool insidePlane(const glm::vec4& v, int plane);
		inline float intersectPlane(const glm::vec4& v1, const glm::vec4& v2, int plane);
		inline ClippedVertex interpolateVertices(const ClippedVertex& v0, const ClippedVertex& v1, float t_Point);
	};
}