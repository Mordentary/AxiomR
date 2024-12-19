#include "renderer.hpp"
#define NOMINMAX
#include <windows.h>
#include <memory>
#include <stdexcept>
#include <vector>
#include <mesh.hpp>
#include <vec.hpp>
#include <stb_image.h>
#include <pipeline.hpp>
#include <shaders/shaders.hpp>

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

	void Renderer::drawMesh(const mat4f& transMat, const Mesh& mesh)
	{
		Vec3f lightDir = { 0, 0, -1 };
		lightDir.normalize();
		m_DefaultPipeline->drawMesh(transMat, mesh);
	}

	void Renderer::drawLine(Vec2f p0, Vec2f p1, Color color)
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
				m_Framebuffer->setPixel(y, x, color);
			}
			else {
				m_Framebuffer->setPixel(x, y, color);
			}
			error2 += derror2;
			if (error2 > dx) {
				y += (y1 > y0 ? 1 : -1);
				error2 -= dx * 2;
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
		m_Framebuffer = std::make_unique<Framebuffer>(m_Width, m_Height, true);
		m_Bitmap = std::make_unique<WindowsBitmap>((HWND)m_WindowHandler, m_Width, m_Height);

		int imageWidth, imageHeight, channels;

		m_Camera = std::make_unique<Camera>(Vec3f{ 0.0,0.0,5.0f }, Vec3f{ 0,0,0 }, Vec3f{ 0,1,0 }, degreeToRad(60), m_Width / m_Height, m_Width, m_Height);
		m_Camera->setViewport(0, 0, m_Width, m_Height);

		m_DefaultPipeline = new Pipeline();
		m_DefaultPipeline->setCamera(m_Camera.get());
		m_DefaultPipeline->setFramebuffer(m_Framebuffer.get());

		FlatShader* flatShader = new FlatShader();
		flatShader->lightDirection = Vec3f(1.0f, -1.0f, 0.5f);
		flatShader->lightDirection.normalize();
		DefaultShader* defaultShader = new DefaultShader();
		defaultShader->lightDirection = Vec3f(0.0f, -1.0f, 0.0f);
		defaultShader->lightDirection.normalize();

		m_DefaultShader = defaultShader;
		m_DefaultPipeline->setShader(m_DefaultShader);
		m_Meshes.emplace_back(std::make_unique<Mesh>("assets/Gas Tank/Gas Tank.obj"));
	}
	Renderer::~Renderer()
	{
		m_Framebuffer.reset();
		m_Bitmap.reset();
		DestroyWindow((HWND)m_WindowHandler);
	}
	void Renderer::render()
	{
		for (auto& mesh : m_Meshes)
		{
			float time = GetTickCount64() / 1000.0f;
			float angleX = 0.5f * time;
			float angleY = 0.0f;
			float angleZ = 0.0f;
			mat4f rotation = mat4f::rotateXYZ(angleX, angleY, angleZ);
			drawMesh(rotation, *mesh);
		}
	}

	void Renderer::run() {
		Color white{ 255,255,255,255 };
		Color red{ 255,0,0,255 };

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
			render();
			// Get client area size
			RECT rect;
			GetClientRect((HWND)m_WindowHandler, &rect);
			int windowWidth = rect.right - rect.left;
			int windowHeight = rect.bottom - rect.top;

			m_Bitmap->copyBuffer(m_Framebuffer->getColorData());
			m_Bitmap->present(windowWidth, windowHeight);

			m_Framebuffer->clearColor({ 0,0,0 });
			m_Framebuffer->clearDepth();

			//Sleep(16);
		}
	}
}