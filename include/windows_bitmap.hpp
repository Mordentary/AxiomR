#pragma once
#define NOMINMAX
#include <windows.h>
#include <cstdint>

namespace AR {
	class WindowsBitmap {
	public:
		WindowsBitmap(HWND hwnd, uint32_t width, uint32_t height);
		~WindowsBitmap();

		// Prevent copying
		WindowsBitmap(const WindowsBitmap&) = delete;
		WindowsBitmap& operator=(const WindowsBitmap&) = delete;

		// Allow moving
		WindowsBitmap(WindowsBitmap&& other) noexcept;
		WindowsBitmap& operator=(WindowsBitmap&& other) noexcept;

		void copyBuffer(const void* data);
		void bitBlit(int windowWidth, int windowHeight, bool stretch);

		// Getters
		void* getBits() const { return m_Bits; }
		uint32_t getWidth() const { return m_Width; }
		uint32_t getHeight() const { return m_Height; }

	private:
		HDC m_DeviceContext;
		HDC m_MemoryDC;
		HBITMAP m_Bitmap;
		HBITMAP m_OldBitmap;
		void* m_Bits;
		uint32_t m_Width;
		uint32_t m_Height;
		HWND m_Window;
	};
}