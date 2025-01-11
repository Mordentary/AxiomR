#include "windows_bitmap.hpp"
#include <stdexcept>
#include <tracy\Tracy.hpp>

namespace AR {
	WindowsBitmap::WindowsBitmap(HWND hwnd, uint32_t width, uint32_t height)
		: m_Width(width)
		, m_Height(height)
		, m_Bits(nullptr)
		, m_Window(hwnd)
	{
		m_DeviceContext = GetDC(hwnd);
		if (!m_DeviceContext) {
			throw std::runtime_error("Failed to get device context");
		}

		BITMAPINFO bmi = {};
		bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
		bmi.bmiHeader.biWidth = width;
		bmi.bmiHeader.biHeight = height;
		bmi.bmiHeader.biPlanes = 1;
		bmi.bmiHeader.biBitCount = 32;
		bmi.bmiHeader.biCompression = BI_RGB;

		m_Bitmap = CreateDIBSection(m_DeviceContext, &bmi, DIB_RGB_COLORS, &m_Bits, nullptr, 0);
		if (!m_Bitmap) {
			ReleaseDC(m_Window, m_DeviceContext);
			throw std::runtime_error("Failed to create DIB section");
		}

		m_MemoryDC = CreateCompatibleDC(m_DeviceContext);
		if (!m_MemoryDC) {
			DeleteObject(m_Bitmap);
			ReleaseDC(m_Window, m_DeviceContext);
			throw std::runtime_error("Failed to create memory DC");
		}

		m_OldBitmap = (HBITMAP)SelectObject(m_MemoryDC, m_Bitmap);
		SetStretchBltMode(m_DeviceContext, HALFTONE);
	}

	WindowsBitmap::~WindowsBitmap() {
		if (m_OldBitmap) SelectObject(m_MemoryDC, m_OldBitmap);
		if (m_Bitmap) DeleteObject(m_Bitmap);
		if (m_MemoryDC) DeleteDC(m_MemoryDC);
		if (m_DeviceContext) ReleaseDC(m_Window, m_DeviceContext);
	}

	WindowsBitmap::WindowsBitmap(WindowsBitmap&& other) noexcept
		: m_DeviceContext(other.m_DeviceContext)
		, m_MemoryDC(other.m_MemoryDC)
		, m_Bitmap(other.m_Bitmap)
		, m_OldBitmap(other.m_OldBitmap)
		, m_Bits(other.m_Bits)
		, m_Width(other.m_Width)
		, m_Height(other.m_Height)
		, m_Window(other.m_Window)
	{
		other.m_DeviceContext = nullptr;
		other.m_MemoryDC = nullptr;
		other.m_Bitmap = nullptr;
		other.m_OldBitmap = nullptr;
		other.m_Bits = nullptr;
		other.m_Window = nullptr;
	}

	WindowsBitmap& WindowsBitmap::operator=(WindowsBitmap&& other) noexcept {
		if (this != &other) {
			// Clean up existing resources
			if (m_OldBitmap) SelectObject(m_MemoryDC, m_OldBitmap);
			if (m_Bitmap) DeleteObject(m_Bitmap);
			if (m_MemoryDC) DeleteDC(m_MemoryDC);
			if (m_DeviceContext) ReleaseDC(m_Window, m_DeviceContext);

			// Move resources
			m_DeviceContext = other.m_DeviceContext;
			m_MemoryDC = other.m_MemoryDC;
			m_Bitmap = other.m_Bitmap;
			m_OldBitmap = other.m_OldBitmap;
			m_Bits = other.m_Bits;
			m_Width = other.m_Width;
			m_Height = other.m_Height;
			m_Window = other.m_Window;

			// Null out other's resources
			other.m_DeviceContext = nullptr;
			other.m_MemoryDC = nullptr;
			other.m_Bitmap = nullptr;
			other.m_OldBitmap = nullptr;
			other.m_Bits = nullptr;
			other.m_Window = nullptr;
		}
		return *this;
	}

	void WindowsBitmap::copyBuffer(const void* data) {
		ZoneScoped;
		memcpy(m_Bits, data, m_Width * m_Height * 4);
	}

	void WindowsBitmap::bitBlit(int windowWidth, int windowHeight, bool stretchNecessary) {
		ZoneScoped;
		if (!stretchNecessary) {
			BitBlt(
				m_DeviceContext,
				0, 0,
				m_Width, m_Height,
				m_MemoryDC,
				0, 0,
				SRCCOPY
			);
		}
		else
		{
			SetStretchBltMode(m_DeviceContext, COLORONCOLOR);

			StretchBlt(
				m_DeviceContext,
				0, 0,
				windowWidth, windowHeight,
				m_MemoryDC,
				0, 0,
				m_Width, m_Height,
				SRCCOPY
			);
		}
	}
}