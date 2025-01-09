#pragma once
#include <string>
#include "math.hpp"

namespace AR {
	class Texture {
	public:
		Texture(const std::string& filePath);
		~Texture();

		bool loadFromFile(const std::string& filePath);
		inline glm::vec4 sample(const glm::vec2& uv) const
		{
			if (!m_Data) {
				return { 0, 0, 0, 255 }; // Default color if texture is not loaded
			}
			int x = static_cast<int>(uv.x * (m_Width - 1));
			int y = static_cast<int>(uv.y * (m_Height - 1));

			x = std::max(0, std::min(x, m_Width - 1));
			y = std::max(0, std::min(y, m_Height - 1));

			y = m_Height - 1 - y;

			return getPixel(x, y);
		}
		inline glm::vec4 getPixel(int x, int y) const {
			if (x < 0 || x >= m_Width || y < 0 || y >= m_Height) {
				return glm::vec4{ 0.0f, 0.0f, 0.0f, 0.0f };
			}
			unsigned char* pixel = m_Data + (y * m_Width + x) * 4;

			float r = 0.0f, g = 0.0f, b = 0.0f, a = 255.0f; // Default values
			if (m_Channels > 0) r = pixel[0];       // Red
			if (m_Channels > 1) g = pixel[1];       // Green
			if (m_Channels > 2) b = pixel[2];       // Blue
			if (m_Channels > 3) a = pixel[3];       // Alpha

			float invColorConst = 1.0f / 255.f;
			return glm::vec4{
				r * invColorConst,
				g * invColorConst,
				b * invColorConst,
				a * invColorConst
			};
		}

		inline int getWidth() const { return m_Width; }
		inline int getHeight() const { return m_Height; }
		inline int getChannels() const { return m_Channels; }

		// Prevent copy construction and assignment
		Texture(const Texture&) = delete;
		Texture& operator=(const Texture&) = delete;

		// Allow move construction and assignment
		Texture(Texture&&) = default;
		Texture& operator=(Texture&&) = default;

	private:
		unsigned char* m_Data = nullptr;
		int m_Width;
		int m_Height;
		int m_Channels;
	};
}