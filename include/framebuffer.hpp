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

		void setPixel(int x, int y, const Color& color);
		Color getPixel(int x, int y) const;

		void setDepth(int x, int y, float depthValue);
		float getDepth(int x, int y) const;

		void clearColor(const Color& color);
		void clearDepth(float depthValue = std::numeric_limits<float>::infinity());

		const unsigned char* getColorData() const;
		const float* getDepthData() const;

		uint32_t getWidth() const;
		uint32_t getHeight() const;
		bool isDepthBufferEnabled() const;

	private:
		bool inBounds(int x, int y) const;

		int m_Width;
		int m_Height;
		bool m_UseDepth;
		std::vector<unsigned char> m_Data;    // Color data: RGBA per pixel
		std::vector<float> m_DepthData;       // Depth data: 1 float per pixel
	};
}