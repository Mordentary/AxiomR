// Pipeline.hpp
#pragma once

#include <algorithm>
#include <limits>
#include "shader_interface.hpp"
#include "camera.hpp"
#include "renderer.hpp"
#include "mesh.hpp"

namespace AR
{
	class Pipeline {
		friend IShader;
	public:
		Pipeline();

		void setShader(IShader* shader);
		void setCamera(const Camera* cam);
		void setFramebuffer(Framebuffer* fb);
		void drawMesh(const mat4f& modelMatrix, const Mesh& mesh);
		mat4f getMVPMat() { return m_MVP; };
		mat4f getViewportMat() { return m_Camera->getViewportMatrix(); }
	private:
		IShader* m_Shader;
		const Camera* m_Camera;
		Framebuffer* m_Framebuffer;
		mat4f m_MVP;
		Vec3f barycentric(const Vec2f& A, const Vec2f& B, const Vec2f& C, const Vec2f& P) const;
		void rasterizeTriangle(Vec4f clip[3]);
	};
}