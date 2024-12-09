#pragma once
#include "vector"
#include"vec.hpp"
#include "camera.hpp"
namespace AR
{
	struct Color {
		unsigned char r, g, b, a;
	};
	class Framebuffer {
	public:
		Framebuffer(int width, int height, bool useDepth = false)
			: m_Width(width), m_Height(height), m_UseDepth(useDepth)
		{
			// Each pixel has 4 components: R, G, B, A
			m_Data.resize(static_cast<size_t>(m_Width) * m_Height * 4, 0);

			if (m_UseDepth) {
				m_DepthData.resize(static_cast<size_t>(m_Width) * m_Height, 1.0f);
			}
		}

		// Enable or disable depth buffer after construction
		void setUseDepthBuffer(bool useDepth) {
			m_UseDepth = useDepth;
			if (m_UseDepth && m_DepthData.empty()) {
				m_DepthData.resize(static_cast<size_t>(m_Width) * m_Height, 1.0f);
			}
			else if (!m_UseDepth) {
				m_DepthData.clear();
			}
		}

		void setPixel(int x, int y, const Color& color) {
			if (!inBounds(x, y)) return;
			int index = (y * m_Width + x) * 4;
			m_Data[index + 0] = color.b;
			m_Data[index + 1] = color.g;
			m_Data[index + 2] = color.r;
			m_Data[index + 3] = color.a;
		}

		Color getPixel(int x, int y) const {
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

		void setDepth(int x, int y, float depthValue) {
			if (!m_UseDepth) return;
			if (!inBounds(x, y)) return;
			m_DepthData[y * m_Width + x] = depthValue;
		}

		float getDepth(int x, int y) const {
			if (!m_UseDepth) return 1.0f;
			if (!inBounds(x, y)) return 1.0f;
			return m_DepthData[y * m_Width + x];
		}

		void clearColor(const Color& color) {
			for (int y = 0; y < m_Height; y++) {
				for (int x = 0; x < m_Width; x++) {
					setPixel(x, y, color);
				}
			}
		}

		void clearDepth(float depthValue = 1.0f) {
			if (!m_UseDepth) return;
			for (size_t i = 0; i < m_DepthData.size(); i++) {
				m_DepthData[i] = depthValue;
			}
		}

		const unsigned char* getColorData() const {
			return m_Data.data();
		}

		const float* getDepthData() const {
			return m_UseDepth ? m_DepthData.data() : nullptr;
		}

		uint32_t getWidth() const { return m_Width; }
		uint32_t getHeight() const { return m_Height; }
		bool isDepthBufferEnabled() const { return m_UseDepth; }

	private:
		bool inBounds(int x, int y) const {
			return x >= 0 && x < m_Width && y >= 0 && y < m_Height;
		}

	private:
		int m_Width;
		int m_Height;
		bool m_UseDepth;
		std::vector<unsigned char> m_Data;    // Color data: RGBA per pixel
		std::vector<float> m_DepthData;       // Depth data: 1 float per pixel
	};
}