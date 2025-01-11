#include "Framebuffer.hpp"
#include <tracy\Tracy.hpp>

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

	void Framebuffer::clearColor(const Color& color) {
		ZoneScoped;

		uint32_t packedColor = (static_cast<uint32_t>(color.a) << 24) |
			(static_cast<uint32_t>(color.r) << 16) |
			(static_cast<uint32_t>(color.g) << 8) |
			static_cast<uint32_t>(color.b);
		size_t totalPixels = static_cast<size_t>(m_Width) * static_cast<size_t>(m_Height);
		uint32_t* pixelData = reinterpret_cast<uint32_t*>(m_Data.data());
		std::fill_n(pixelData, totalPixels, packedColor);
	}

	void Framebuffer::clearDepth(float depthValue) {
		ZoneScoped;
		if (!m_UseDepth) return;
		std::fill(m_DepthData.begin(), m_DepthData.end(), depthValue);
	}

	const uint8_t* Framebuffer::getColorData() const {
		return m_Data.data();
	}

	const float* Framebuffer::getDepthData() const {
		return m_UseDepth ? m_DepthData.data() : nullptr;
	}

	uint8_t* Framebuffer::getColorData() {
		return m_Data.data();
	}

	float* Framebuffer::getDepthData() {
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
}