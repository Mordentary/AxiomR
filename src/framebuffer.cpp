#include "Framebuffer.hpp"

namespace AR
{
	Framebuffer::Framebuffer(int width, int height, bool useDepth)
		: m_Width(width), m_Height(height), m_UseDepth(useDepth)
	{
		m_Data.resize(static_cast<size_t>(m_Width) * m_Height * 4, 0);

		if (m_UseDepth) {
			m_DepthData.resize(static_cast<size_t>(m_Width) * m_Height, std::numeric_limits<float>::infinity());
		}
	}

	void Framebuffer::setUseDepthBuffer(bool useDepth) {
		m_UseDepth = useDepth;
		if (m_UseDepth && m_DepthData.empty()) {
			m_DepthData.resize(static_cast<size_t>(m_Width) * m_Height, 1.0f);
		}
		else if (!m_UseDepth) {
			m_DepthData.clear();
		}
	}

	void Framebuffer::setPixel(int x, int y, const Color& color) {
		if (!inBounds(x, y)) return;
		int index = (y * m_Width + x) * 4;
		m_Data[index + 0] = color.b;
		m_Data[index + 1] = color.g;
		m_Data[index + 2] = color.r;
		m_Data[index + 3] = color.a;
	}

	Color Framebuffer::getPixel(int x, int y) const {
		if (!inBounds(x, y)) {
			return { 0, 0, 0, 0 };
		}
		int index = (y * m_Width + x) * 4;
		return {
			m_Data[index + 2], // R
			m_Data[index + 1], // G
			m_Data[index + 0], // B
			m_Data[index + 3]  // A
		};
	}

	void Framebuffer::setDepth(int x, int y, float depthValue) {
		if (!m_UseDepth) return;
		if (!inBounds(x, y)) return;
		m_DepthData[y * m_Width + x] = depthValue;
	}

	float Framebuffer::getDepth(int x, int y) const {
		if (!m_UseDepth) return std::numeric_limits<float>::infinity();
		if (!inBounds(x, y)) return std::numeric_limits<float>::infinity();
		return m_DepthData[y * m_Width + x];
	}

	void Framebuffer::clearColor(const Color& color) {
		for (int y = 0; y < m_Height; y++) {
			for (int x = 0; x < m_Width; x++) {
				setPixel(x, y, color);
			}
		}
	}

	void Framebuffer::clearDepth(float depthValue) {
		if (!m_UseDepth) return;
		std::fill(m_DepthData.begin(), m_DepthData.end(), depthValue);
	}

	const unsigned char* Framebuffer::getColorData() const {
		return m_Data.data();
	}

	const float* Framebuffer::getDepthData() const {
		return m_UseDepth ? m_DepthData.data() : nullptr;
	}

	uint32_t Framebuffer::getWidth() const {
		return static_cast<uint32_t>(m_Width);
	}

	uint32_t Framebuffer::getHeight() const {
		return static_cast<uint32_t>(m_Height);
	}

	bool Framebuffer::isDepthBufferEnabled() const {
		return m_UseDepth;
	}

	bool Framebuffer::inBounds(int x, int y) const {
		return x >= 0 && x < m_Width && y >= 0 && y < m_Height;
	}
}