#pragma once
#include <vector>
#include "windows_bitmap.hpp"

#include"framebuffer.hpp"
#include"vec.hpp"
#include "camera.hpp"

namespace AR
{
	class Vertex;
	class Mesh;

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
		std::unique_ptr<Framebuffer> m_Framebuffer;
		uint32_t m_Width, m_Height;
		void* m_WindowHandler;
		std::unique_ptr<WindowsBitmap> m_Bitmap;
		std::unique_ptr<Camera> m_Camera;

		int m_ImageWidth = 0;
		int m_ImageHeight = 0;
		unsigned char* m_ImageData = nullptr;
	};
}