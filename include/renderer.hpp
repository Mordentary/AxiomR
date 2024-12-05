#pragma once
#include "vector"
#include "vec.hpp"
#include "windows_bitmap.hpp"
namespace AR
{
	class Vertex;
	class Mesh;
	struct Color {
		unsigned char r, g, b, a;
	};

	struct Point3 {
		float x = 0, y = 0, z = 0;
	};

	struct Point2 {
		float x = 0, y = 0;
	};

	struct Buffer {
		Buffer(int w, int h) : width(w), height(h) {
			data.resize(w * h * 4, 0);
		}

		void setPixel(int x, int y, Color color) {
			if (x < 0 || x >= width || y < 0 || y >= height) return;  // Bounds checking
			int index = (y * width + x) * 4;
			data[index + 0] = color.b;    // Blue
			data[index + 1] = color.g;    // Green
			data[index + 2] = color.r;    // Red
			data[index + 3] = color.a;    // Alpha
		}

		Color getPixel(int x, int y) const {
			if (x < 0 || x >= width || y < 0 || y >= height)
				return { 0, 0, 0, 0 };  // Return transparent black for out of bounds

			int actualY = height - y - 1;
			int index = (actualY * width + x) * 4;

			return {
				data[index + 2],  // Red
				data[index + 1],  // Green
				data[index + 0],  // Blue
				data[index + 3]   // Alpha
			};
		}

		// Get raw data pointer for rendering
		const unsigned char* getData() const {
			return data.data();
		}

		uint32_t getWidth() const {
			return width;
		}

		uint32_t getHeight() const {
			return height;
		}

		// Clear buffer to a specific color
		void clear(Color color) {
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					setPixel(x, y, color);
				}
			}
		}

	private:
		std::vector<unsigned char> data;
		int width;
		int height;
	};

	class Renderer
	{
	public:
		explicit Renderer() = default;
		~Renderer();
		void init(uint32_t screenWidth, uint32_t screenHeight);
		void run();
		void drawLine(Point2 p1, Point2 p2, Color color);
		void drawTriangle(const Vertex& p0, const  Vertex& p1, const  Vertex& p2, Vec3f normal, Vec3f ligthDir);
		void drawMesh(const Mesh& mesh);
	private:
		std::unique_ptr<Buffer> m_ScreenBuffer;
		uint32_t m_Width, m_Height;
		void* m_WindowHandler;
		std::unique_ptr<WindowsBitmap> m_Bitmap;
		std::vector<float> m_DepthBuffer;

		int m_ImageWidth = 0;
		int m_ImageHeight= 0;
		unsigned char* m_ImageData = nullptr;
	};
}