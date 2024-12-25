#pragma once
#include <string>
#include "math.hpp"

namespace AR {
	class Texture {
	public:
		Texture(const std::string& filePath);
		~Texture();

		bool loadFromFile(const std::string& filePath);
		Vec4f sample(const Vec2f& uv) const;
		Vec4f getPixel(int x, int y) const;

		int getWidth() const { return m_Width; }
		int getHeight() const { return m_Height; }
		int getChannels() const { return m_Channels; }

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