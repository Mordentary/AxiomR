#include "Renderer.hpp"

#define NOMINMAX
#include <windows.h>
#include <memory>
#include <stdexcept>
#include <vector>
#include <mesh.hpp>
#include <vec.hpp>

namespace AR {
	const char CLASS_NAME[] = "AxiomR Window";
	const char TITLE[] = "AxiomR";

	LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
		switch (uMsg) {
		case WM_DESTROY:
			PostQuitMessage(0);
			return 0;
		}
		return DefWindowProc(hwnd, uMsg, wParam, lParam);
	}

	void Renderer::drawMesh(const Mesh& mesh)
	{
		Vec3f lightDir = { 0,0,-1 };
		auto& meshVertices = mesh.getVertices();
		for (auto& face : mesh.getFaces())
		{
			uint32_t faceSize = face.vertexIndices.size();
			for (int i = 0; i < faceSize; i += 3)
			{
				const Vertex& v1 = meshVertices[face.vertexIndices[i]];
				const Vertex& v2 = meshVertices[face.vertexIndices[i + 1]];
				const Vertex& v3 = meshVertices[face.vertexIndices[i + 2]];
				Vec3f normal = ((v3.position - v1.position).cross(v2.position - v1.position));
				normal.normalize();
				lightDir.normalize();
				drawTriangle(v1, v2, v3, normal, lightDir);
			}
		}
	}

	void Renderer::drawLine(Point2 p0, Point2 p1, Color color)
	{
		int x0 = p0.x, y0 = p0.y;
		int x1 = p1.x, y1 = p1.y;

		bool steep = false;
		if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
			std::swap(x0, y0);
			std::swap(x1, y1);
			steep = true;
		}
		if (x0 > x1) {
			std::swap(x0, x1);
			std::swap(y0, y1);
		}

		int dx = x1 - x0;
		int dy = y1 - y0;
		int derror2 = std::abs(dy) * 2;
		int error2 = 0;
		int y = y0;

		for (int x = x0; x <= x1; x++) {
			if (steep) {
				m_ScreenBuffer->setPixel(y, x, color);
			}
			else {
				m_ScreenBuffer->setPixel(x, y, color);
			}
			error2 += derror2;
			if (error2 > dx) {
				y += (y1 > y0 ? 1 : -1);
				error2 -= dx * 2;
			}
		}
	}

	void Renderer::drawTriangle(const Vertex& p0, const  Vertex& p1, const  Vertex& p2, Vec3f normal, Vec3f lightDir)
	{
		// Get viewport dimensions
		int width = m_ScreenBuffer->getWidth();
		int height = m_ScreenBuffer->getHeight();

		// Transform from NDC [-1,1] to screen space [0,width/height]
		Vec2i pts[3] = {
			{static_cast<int>((p0.position.x + 1.0f) * width * 0.5f),
			 static_cast<int>((p0.position.y + 1.0f) * height * 0.5f)},
			{static_cast<int>((p1.position.x + 1.0f) * width * 0.5f),
			 static_cast<int>((p1.position.y + 1.0f) * height * 0.5f)},
			{static_cast<int>((p2.position.x + 1.0f) * width * 0.5f),
			 static_cast<int>((p2.position.y + 1.0f) * height * 0.5f)}
		};

		// Calculate bounding box
		Vec2i bboxmin(m_ScreenBuffer->getWidth() - 1, m_ScreenBuffer->getHeight() - 1);
		Vec2i bboxmax(0, 0);
		Vec2i clamp(m_ScreenBuffer->getWidth() - 1, m_ScreenBuffer->getHeight() - 1);

		// Find min/max points for bounding box
		for (int i = 0; i < 3; i++) {
			bboxmin.x = std::max(0, std::min(bboxmin.x, pts[i].x));
			bboxmin.y = std::max(0, std::min(bboxmin.y, pts[i].y));
			bboxmax.x = std::min(clamp.x, std::max(bboxmax.x, pts[i].x));
			bboxmax.y = std::min(clamp.y, std::max(bboxmax.y, pts[i].y));
		}

		// Rasterize
		Vec2i P;
		for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
			for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
				Vec3f bc_screen = barycentric(pts, P);
				if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;

				// Interpolate Z-coordinate for depth
				float z = p0.position.z * bc_screen.x + p1.position.z * bc_screen.y + p2.position.z * bc_screen.z;

				// Interpolate normals
				//Vec3f normal = {
				//	p0.normal.x * bc_screen.x + p1.normal.x * bc_screen.y + p2.normal.x * bc_screen.z,
				//	p0.normal.y * bc_screen.x + p1.normal.y * bc_screen.y + p2.normal.y * bc_screen.z,
				//	p0.normal.z * bc_screen.x + p1.normal.z * bc_screen.y + p2.normal.z * bc_screen.z
				//};

				// Calculate lighting intensity
				float intensity = normal.dot(lightDir);

				if (intensity > 0)
				{
					Color color = {
						static_cast<uint8_t>(255 * intensity),
						static_cast<uint8_t>(255 * intensity),
						static_cast<uint8_t>(255 * intensity),
						255
					};

					//if (z < m_DepthBuffer[P.x + P.y * width]) {
					//	m_DepthBuffer[P.x + P.y * width] = z;
					//}
					m_ScreenBuffer->setPixel(P.x, P.y, color);
				}
			}
		}
	}

	void Renderer::init(uint32_t screenWidth, uint32_t screenHeight)
	{
		m_Width = screenWidth;
		m_Height = screenHeight;

		WNDCLASS wc = {};
		wc.lpfnWndProc = WindowProc;
		wc.hInstance = GetModuleHandle(nullptr);
		wc.lpszClassName = CLASS_NAME;
		RegisterClass(&wc);

		// Create window
		HWND hwnd = CreateWindowEx(
			0,
			CLASS_NAME,
			TITLE,
			WS_OVERLAPPEDWINDOW,
			CW_USEDEFAULT, CW_USEDEFAULT,
			m_Width, m_Height,
			nullptr,
			nullptr,
			GetModuleHandle(nullptr),
			nullptr
		);

		if (!hwnd) {
			throw std::runtime_error("Failed to create window");
		}
		m_WindowHandler = hwnd;
		// Show window
		ShowWindow(hwnd, SW_SHOW);

		//TODO: Loading resources not here
		Mesh mesh("assets/african_head.obj");
		m_ScreenBuffer = std::make_unique<Buffer>(m_Width, m_Height);
		m_Bitmap = std::make_unique<WindowsBitmap>((HWND)m_WindowHandler, m_Width, m_Height);
		m_DepthBuffer.resize(m_Width * m_Height, std::numeric_limits<float>::infinity());
		drawMesh(mesh);
	}
	Renderer::~Renderer()
	{
		m_ScreenBuffer.reset();
		m_Bitmap.reset();
		DestroyWindow((HWND)m_WindowHandler);
	}

	void Renderer::run() {
		Color white{ 255,255,255,255 };
		Color red{ 255,0,0,255 };

		// Message loop
		bool running = true;
		while (running) {
			MSG msg = {};
			while (PeekMessage(&msg, nullptr, 0, 0, PM_REMOVE)) {
				if (msg.message == WM_QUIT) {
					running = false;
					break;
				}
				TranslateMessage(&msg);
				DispatchMessage(&msg);
			}

			// Get client area size
			RECT rect;
			GetClientRect((HWND)m_WindowHandler, &rect);
			int windowWidth = rect.right - rect.left;
			int windowHeight = rect.bottom - rect.top;

			m_Bitmap->copyBuffer(m_ScreenBuffer->getData());
			m_Bitmap->render(windowWidth, windowHeight);

			Sleep(16);
		}
	}
}