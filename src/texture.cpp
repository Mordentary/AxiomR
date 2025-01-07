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
}