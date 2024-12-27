#pragma once
#define NOMINMAX
#include <windows.h>
#include <memory>
#include "framebuffer.hpp"
#include "windows_bitmap.hpp"

namespace AR {
	class Renderer;

	class Window {
		friend class Renderer;
	public:
		Window(uint32_t width, uint32_t height, const char* title);
		~Window();

		HWND getHandle() const { return m_WindowHandle; }
		uint32_t getWidth() const { return m_Width; }
		uint32_t getHeight() const { return m_Height; }
		glm::vec2 getMousePos() const { return m_MousePos; }
		bool isMouseButtonDown(int button) const; // 0 for left, 1 for right, 2 for middle
		bool isKeyDown(char key) const;

		void show();
		void processMessages();
		bool isRunning() const { return m_Running; }

		void present(const Framebuffer& framebuffer);

	private:
		static LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam);
		// Internal helper to handle messages for this specific instance
		LRESULT handleMessage(UINT uMsg, WPARAM wParam, LPARAM lParam);

	private:
		uint32_t m_Width;
		uint32_t m_Height;
		HWND m_WindowHandle;
		std::unique_ptr<WindowsBitmap> m_Bitmap;
		bool m_Running;
		glm::vec2 m_MousePos;
		bool m_MouseButtons[3];
		bool m_Keys[256];
	};
}
