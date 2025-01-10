#pragma once

#include <vector>
#include <limits>
#include "math.hpp"

namespace AR
{
	class Framebuffer {
	public:
		Framebuffer(int width, int height, bool useDepth = false);
		void setUseDepthBuffer(bool useDepth);

		inline void setPixel(int x, int y, const Color& color) {
			if (!inBounds(x, y)) return;
			int index = (y * m_Width + x) * 4;
			m_Data[index + 0] = color.b;
			m_Data[index + 1] = color.g;
			m_Data[index + 2] = color.r;
			m_Data[index + 3] = color.a;
		}

		inline Color getPixel(int x, int y) const {
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

		inline void setDepth(int x, int y, float depthValue) {
			if (!m_UseDepth) return;
			if (!inBounds(x, y)) return;
			m_DepthData[y * m_Width + x] = depthValue;
		}
		inline float getDepth(int x, int y) const {
			if (!m_UseDepth) return std::numeric_limits<float>::infinity();
			if (!inBounds(x, y)) return std::numeric_limits<float>::infinity();
			return m_DepthData[y * m_Width + x];
		}

		void clearColor(const Color& color);
		void clearDepth(float depthValue = std::numeric_limits<float>::infinity());

		const unsigned char* getColorData() const;
		const float* getDepthData() const;

		uint32_t getWidth() const;
		uint32_t getHeight() const;
		bool isDepthBufferEnabled() const;

	private:
		inline bool inBounds(int x, int y) const {
			return x >= 0 && x < m_Width && y >= 0 && y < m_Height;
		}

		int m_Width;
		int m_Height;
		bool m_UseDepth;
		alignas(32) std::vector<unsigned char> m_Data;    // Color data: RGBA per pixel
		std::vector<float> m_DepthData;       // Depth data: 1 float per pixel
	};
}