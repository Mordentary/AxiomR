#include "texture.hpp"
#include "stb_image.h"
#include <cstdio>

namespace AR {
	Texture::Texture(const std::string& filePath) :
		m_Data(nullptr),
		m_Width(0),
		m_Height(0),
		m_Channels(0)
	{
		loadFromFile(filePath);
	}

	Texture::~Texture() {
		if (m_Data) {
			stbi_image_free(m_Data);
		}
	}

	bool Texture::loadFromFile(const std::string& filePath) {
		int width, height, channels;

		m_Data = stbi_load(filePath.c_str(), &width, &height, &channels, STBI_rgb_alpha);

		if (!m_Data) {
			fprintf(stderr, "Failed to load texture image: %s\n", filePath.c_str());
			return false;
		}

		m_Width = width;
		m_Height = height;
		m_Channels = channels;

		return true;
	}

	Vec4f Texture::getPixel(int x, int y) const {
		// Bounds check
		if (x < 0 || x >= m_Width || y < 0 || y >= m_Height) {
			return Vec4f{ 0.0f, 0.0f, 0.0f, 0.0f }; // Return black with alpha 0 for invalid access
		}

		unsigned char* pixel = m_Data + (y * m_Width + x) * 4;

		float r = 0.0f, g = 0.0f, b = 0.0f, a = 255.0f; // Default values
		if (m_Channels > 0) r = pixel[0];       // Red
		if (m_Channels > 1) g = pixel[1];       // Green
		if (m_Channels > 2) b = pixel[2];       // Blue
		if (m_Channels > 3) a = pixel[3];       // Alpha

		return Vec4f{
			r / 255.0f, // Normalize to [0.0f, 1.0f]
			g / 255.0f,
			b / 255.0f,
			a / 255.0f
		};
	}

	Vec4f Texture::sample(const Vec2f& uv) const
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
}