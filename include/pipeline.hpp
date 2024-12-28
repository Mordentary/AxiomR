#pragma once
#include <algorithm>
#include <limits>
#include "mesh.hpp"
#include <mutex>
namespace AR
{
	class IShader;
	class Camera;
	class Framebuffer;

	class Pipeline {
		friend IShader;
	public:
		Pipeline();
		void setShader(IShader* shader);
		void setCamera(Camera* cam);
		void setFramebuffer(Framebuffer* fb);
		void drawMesh(const glm::mat4& modelMatrix, const Mesh& mesh);
		glm::mat4 getViewportMat();
	private:
		IShader* m_Shader = nullptr;
		Framebuffer* m_Framebuffer = nullptr;
		const Camera* m_Camera = nullptr;
		glm::vec3 barycentric(const glm::vec2& A, const glm::vec2& B, const glm::vec2& C, const glm::vec2& P) const;
		void rasterizeTriangle(const glm::vec4 clip[3]);
		//std::vector<std::array<std::pair<Vertex, glm::vec4>, 3>> clipTriangle(const std::array<std::pair<Vertex, glm::vec4>, 3>& tri);
		void clipTriangle(std::vector<std::pair<Vertex, glm::vec4>>& clippedVertices);
		//std::vector<std::pair<Vertex, glm::vec4>> clipAgainstPlane(const std::vector<std::pair<Vertex, glm::vec4>>& poly, int plane);
		void clipAgainstPlane(std::vector<std::pair<Vertex, glm::vec4>>& poly, int plane);
		inline bool insidePlane(const glm::vec4& v, int plane);
		inline float intersectPlane(const glm::vec4& v1, const glm::vec4& v2, int plane);
		inline std::pair<Vertex, glm::vec4> interpolateVertices(std::pair<Vertex, glm::vec4> v0, std::pair<Vertex, glm::vec4> v1, float t_Point);
	};
}