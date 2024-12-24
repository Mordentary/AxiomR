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
		void setCamera(const Camera* cam);
		void setFramebuffer(Framebuffer* fb);
		void drawMesh(const mat4f& modelMatrix, const Mesh& mesh);
		mat4f getViewportMat();
	private:
		IShader* m_Shader;
		Framebuffer* m_Framebuffer;
		const Camera* m_Camera;
		Vec3f barycentric(const Vec2f& A, const Vec2f& B, const Vec2f& C, const Vec2f& P) const;
		void rasterizeTriangle(const Vec4f clip[3]);
		void clipAndRasterizeTriangle(const Vec4f clip[3]);
		std::vector<std::array<Vec4f, 3>> clipTriangle(const std::array<Vec4f, 3>& tri);
		std::vector<Vec4f> clipAgainstPlane(const std::vector<Vec4f>& poly, int plane);
		bool insidePlane(const Vec4f& v, int plane);
		Vec4f intersectPlane(const Vec4f& v1, const Vec4f& v2, int plane);
	};
}