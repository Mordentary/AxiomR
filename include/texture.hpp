#pragma once
#include <string>
#include "math.hpp"
#include"tracy\Tracy.hpp"
namespace AR {
	class Texture {
	public:
		Texture(const std::string& filePath);
		~Texture();

		bool loadFromFile(const std::string& filePath);
		inline glm::vec4 sample(const glm::vec2& uv) const
		{
			//ZoneScoped;
			if (!m_Data) {
				return { 0, 0, 0, 1 };
			}
			int x = static_cast<int>(uv.x * (m_Width - 1));
			int y = static_cast<int>(uv.y * (m_Height - 1));

			x = std::max(0, std::min(x, m_Width - 1));
			y = std::max(0, std::min(y, m_Height - 1));

			y = m_Height - 1 - y;

			unsigned char* pixel = m_Data + (y * m_Width + x) * 4;
			float inv255 = 1.0f / 255.0f;
			return glm::vec4{
				pixel[0] * inv255,
				pixel[1] * inv255,
				pixel[2] * inv255,
				pixel[3] * inv255
			};
		}
		inline glm::vec4 getPixelRGBA(int x, int y) const {
			unsigned char* pixel = m_Data + (y * m_Width + x) * 4;

			float inv255 = 1.0f / 255.0f;
			return glm::vec4{
				pixel[0] * inv255,
				pixel[1] * inv255,
				pixel[2] * inv255,
				pixel[3] * inv255
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