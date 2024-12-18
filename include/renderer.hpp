#pragma once
#include <vector>
#include "camera.hpp"
#include"framebuffer.hpp"
#include"vec.hpp"
#include "windows_bitmap.hpp"
#include"mesh.hpp"
namespace AR
{
	class Vertex;
	class IShader;
	class Pipeline;

	class Renderer
	{
	public:
		explicit Renderer() = default;
		~Renderer();
		void init(uint32_t screenWidth, uint32_t screenHeight);
		void run();

		void drawLine(Vec2f p1, Vec2f p2, Color color);
		void drawTriangle(const Vertex& p0, const  Vertex& p1, const  Vertex& p2, Vec3f normal, Vec3f ligthDir);
		void drawMesh(const mat4f& trnsfrm, const Mesh& mesh);
	private:
		void render();
	private:
		std::unique_ptr<Framebuffer> m_Framebuffer;
		uint32_t m_Width, m_Height;
		void* m_WindowHandler;
		std::unique_ptr<WindowsBitmap> m_Bitmap;
		std::unique_ptr<Camera> m_Camera;

		Pipeline* m_DefaultPipeline;
		IShader* m_DefaultShader;
		std::vector<std::unique_ptr<Mesh>> m_Meshes;
	};
}