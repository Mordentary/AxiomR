#include "Renderer.hpp"

#define NOMINMAX
#include <windows.h>
#include <memory>
#include <stdexcept>
#include <vector>
#include <mesh.hpp>
#include <vec.hpp>
#include <stb_image.h>

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

	void Renderer::drawMesh(const Mat4f& transMat, const Mesh& mesh)
	{
		Vec3f lightDir = { 0, 0, -1 };
		lightDir.normalize();
		// Get camera matrices
		Mat4f view = m_Camera->getViewMatrix();
		Mat4f proj = m_Camera->getProjectionMatrix();
		Mat4f vp = proj * view;

		auto& vertices = mesh.getVertices();
		for (const auto& face : mesh.getFaces()) {
			for (size_t i = 0; i < face.vertexIndices.size(); i += 3) {
				const Vertex& originalV0 = vertices[face.vertexIndices[i]];
				const Vertex& originalV1 = vertices[face.vertexIndices[i + 1]];
				const Vertex& originalV2 = vertices[face.vertexIndices[i + 2]];

				// Apply rotation to the vertices in real time
				Vec4f hv0 = transMat * Vec4f{ originalV0.position.x, originalV0.position.y, originalV0.position.z, 1.0f };
				Vec4f hv1 = transMat * Vec4f{ originalV1.position.x, originalV1.position.y, originalV1.position.z, 1.0f };
				Vec4f hv2 = transMat * Vec4f{ originalV2.position.x, originalV2.position.y, originalV2.position.z, 1.0f };

				hv0 = vp * hv0;
				hv1 = vp * hv1;
				hv2 = vp * hv2;

				// Perspective divide
				if (hv0.w != 0) { hv0.x /= hv0.w; hv0.y /= hv0.w; hv0.z /= hv0.w; }
				if (hv1.w != 0) { hv1.x /= hv1.w; hv1.y /= hv1.w; hv1.z /= hv1.w; }
				if (hv2.w != 0) { hv2.x /= hv2.w; hv2.y /= hv2.w; hv2.z /= hv2.w; }

				Vertex tv0 = originalV0; tv0.position = { hv0.x, hv0.y, hv0.z };
				Vertex tv1 = originalV1; tv1.position = { hv1.x, hv1.y, hv1.z };
				Vertex tv2 = originalV2; tv2.position = { hv2.x, hv2.y, hv2.z };

				// Calculate face normal in world space if needed
				Vec3f normal = ((tv2.position - tv0.position).cross(tv1.position - tv0.position));
				normal.normalize();

				drawTriangle(tv0, tv1, tv2, normal, lightDir);
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

	void Renderer::drawTriangle(const Vertex& p0, const  Vertex& p1, const  Vertex& p2, Vec3f nor, Vec3f lightDir)
	{
		// Get viewport dimensions
		int width = m_Framebuffer->getWidth();
		int height = m_Framebuffer->getHeight();

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
		Vec2i bboxmin(m_Framebuffer->getWidth() - 1, m_Framebuffer->getHeight() - 1);
		Vec2i bboxmax(0, 0);
		Vec2i clamp(m_Framebuffer->getWidth() - 1, m_Framebuffer->getHeight() - 1);

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
				Vec3f normal = {
					p0.normal.x * bc_screen.x + p1.normal.x * bc_screen.y + p2.normal.x * bc_screen.z,
					p0.normal.y * bc_screen.x + p1.normal.y * bc_screen.y + p2.normal.y * bc_screen.z,
					p0.normal.z * bc_screen.x + p1.normal.z * bc_screen.y + p2.normal.z * bc_screen.z
				};
				normal.normalize();
				normal = -normal;

				Vec2f uv =
				{
					(p0.uv.x * bc_screen.x + p1.uv.x * bc_screen.y + p2.uv.x * bc_screen.z),
					(p0.uv.y * bc_screen.x + p1.uv.y * bc_screen.y + p2.uv.y * bc_screen.z),
				};
				int texX = static_cast<int>(uv.x * (m_ImageWidth - 1));
				int texY = static_cast<int>(uv.y * (m_ImageHeight - 1));
				texY = m_ImageHeight - texY - 1;
				texX = std::max(0, std::min(texX, m_ImageWidth - 1));
				texY = std::max(0, std::min(texY, m_ImageHeight - 1));

				// Calculate lighting intensity
				float intensity = normal.dot(lightDir);

				if (intensity > 0)
				{
					int pixelIndex = (texY * m_ImageWidth + texX) * 3;
					unsigned char red = m_ImageData[pixelIndex + 0];
					unsigned char green = m_ImageData[pixelIndex + 1];
					unsigned char blue = m_ImageData[pixelIndex + 2];
					Color color = {
						static_cast<uint8_t>(red * intensity),
						static_cast<uint8_t>(green * intensity),
						static_cast<uint8_t>(blue * intensity),
						255
					};

					if (z > m_Framebuffer->getDepth(P.x, P.y)) {
						m_Framebuffer->setDepth(P.x, P.y, z);
						m_Framebuffer->setPixel(P.x, P.y, color);
					}
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
		m_Framebuffer = std::make_unique<Framebuffer>(m_Width, m_Height, true);
		m_Bitmap = std::make_unique<WindowsBitmap>((HWND)m_WindowHandler, m_Width, m_Height);

		int imageWidth, imageHeight, channels;

		m_ImageData = stbi_load("assets/african_head_diffuse.tga", &imageWidth, &imageHeight, &channels, 0);

		if (!m_ImageData)
			fprintf(stderr, "Failed to load texture image!\n");
		m_ImageWidth = imageWidth;
		m_ImageHeight = imageHeight;
		m_Camera = std::make_unique<Camera>();
	}
	Renderer::~Renderer()
	{
		m_Framebuffer.reset();
		m_Bitmap.reset();
		DestroyWindow((HWND)m_WindowHandler);
	}

	void Renderer::run() {
		Color white{ 255,255,255,255 };
		Color red{ 255,0,0,255 };

		Mesh mesh("assets/african_head.obj");

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

			float time = GetTickCount64() / 1000.0f;
			float angleX = 0.2f;
			float angleY = 0.6f * (time);
			float angleZ = 0.0f;
			Mat4f rotation = Mat4f::rotateXYZ(angleX, angleY, angleZ);
			drawMesh(rotation, mesh);

			// Get client area size
			RECT rect;
			GetClientRect((HWND)m_WindowHandler, &rect);
			int windowWidth = rect.right - rect.left;
			int windowHeight = rect.bottom - rect.top;

			m_Bitmap->copyBuffer(m_Framebuffer->getColorData());
			m_Bitmap->render(windowWidth, windowHeight);

			m_Framebuffer->clearColor({ 0,0,0 });
			m_Framebuffer->clearDepth();

			//Sleep(16);
		}
	}
}