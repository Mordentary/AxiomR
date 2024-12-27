#include "window.hpp"
#include <windowsx.h>
#include <tracy\Tracy.hpp>

namespace AR {
	const char CLASS_NAME[] = "AxiomR Window Class";

	Window::Window(uint32_t width, uint32_t height, const char* title)
		: m_Width(width), m_Height(height), m_WindowHandle(nullptr), m_Running(true)
	{
		for (int i = 0; i < 3; ++i) {
			m_MouseButtons[i] = false;
		}

		for (int i = 0; i < 256; ++i) {
			m_Keys[i] = false;
		}

		WNDCLASS wc = {};
		wc.lpfnWndProc = WindowProc;
		wc.hInstance = GetModuleHandle(nullptr);
		wc.lpszClassName = CLASS_NAME;
		RegisterClass(&wc);

		// Create the window
		m_WindowHandle = CreateWindowEx(
			0,
			CLASS_NAME,
			title,
			WS_OVERLAPPEDWINDOW,
			CW_USEDEFAULT, CW_USEDEFAULT,
			m_Width, m_Height,
			nullptr,
			nullptr,
			GetModuleHandle(nullptr),
			this
		);

		if (!m_WindowHandle) {
			throw std::runtime_error("Failed to create window");
		}

		SetWindowLongPtr(m_WindowHandle, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(this));

		// Initialize the bitmap
		m_Bitmap = std::make_unique<WindowsBitmap>(m_WindowHandle, m_Width, m_Height);
	}

	Window::~Window() {
		m_Bitmap.reset();
		DestroyWindow(m_WindowHandle);
	}

	void Window::show() {
		ShowWindow(m_WindowHandle, SW_SHOW);
	}

	void Window::processMessages() {
		MSG msg = {};
		while (PeekMessage(&msg, nullptr, 0, 0, PM_REMOVE)) {
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
	}

	void Window::present(const Framebuffer& framebuffer)
	{

		ZoneScoped;
		// Get client area size
		RECT rect;
		GetClientRect(m_WindowHandle, &rect);
		int windowWidth = rect.right - rect.left;
		int windowHeight = rect.bottom - rect.top;

		m_Bitmap->copyBuffer(framebuffer.getColorData());
		m_Bitmap->present(windowWidth, windowHeight);
	}

	bool Window::isKeyDown(char key) const {
		return m_Keys[static_cast<unsigned char>(key)];
	}
	LRESULT CALLBACK Window::WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
		Window* pWindow = reinterpret_cast<Window*>(GetWindowLongPtr(hwnd, GWLP_USERDATA));

		if (pWindow) {
			return pWindow->handleMessage(uMsg, wParam, lParam);
		}
		else {
			return DefWindowProc(hwnd, uMsg, wParam, lParam);
		}
	}

	bool Window::isMouseButtonDown(int button) const {
		if (button >= 0 && button < 3) {
			return m_MouseButtons[button];
		}
		return false;
	}

	LRESULT Window::handleMessage(UINT uMsg, WPARAM wParam, LPARAM lParam) {
		switch (uMsg) {
		case WM_DESTROY:
			m_Running = false;
			PostQuitMessage(0);
			return 0;
		case WM_MOUSEMOVE: {
			m_MousePos.x = (float)GET_X_LPARAM(lParam);
			m_MousePos.y = (float) -GET_Y_LPARAM(lParam);
			return 0;
		}
		case WM_LBUTTONDOWN: {
			m_MouseButtons[0] = true;
			SetCapture(m_WindowHandle);
			return 0;
		}
		case WM_LBUTTONUP: {
			m_MouseButtons[0] = false;
			ReleaseCapture();
			return 0;
		}
		case WM_RBUTTONDOWN: {
			m_MouseButtons[1] = true;
			return 0;
		}
		case WM_RBUTTONUP: {
			m_MouseButtons[1] = false;
			return 0;
		}
		case WM_MBUTTONDOWN: {
			m_MouseButtons[2] = true;
			return 0;
		}
		case WM_MBUTTONUP: {
			m_MouseButtons[2] = false;
			return 0;
		}
		case WM_KEYDOWN: {
			m_Keys[wParam] = true;
			return 0;
		}
		case WM_KEYUP: {
			m_Keys[wParam] = false;
			return 0;
		}
		case WM_SIZE: {
			//todo:resize
			//m_Width = LOWORD(lParam);
			//m_Height = HIWORD(lParam);

			//if (m_Renderer) {
			//	//m_Renderer->resize(m_Width, m_Height);
			//}
			return 0;
		}
		default:
			return DefWindowProc(m_WindowHandle, uMsg, wParam, lParam);
		}
	}
}
