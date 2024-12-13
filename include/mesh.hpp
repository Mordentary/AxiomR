#pragma once
#include <vector>
#include <unordered_map>
#include"vec.hpp"

namespace AR {
	struct Vertex {
		Vec3f position;
		Vec2f uv;
		Vec3f normal;
		bool operator==(const Vertex& other) const {
			return position == other.position && uv == other.uv && normal == other.normal;
		}
	};

	struct VertexHash {
		std::size_t operator()(const Vertex& v) const {
			// Simple hashing approach:
			std::hash<float> hf;
			size_t h = 0;
			auto hash_combine = [&](size_t& seed, float val) {
				size_t hVal = hf(val);
				seed ^= hVal + 0x9e3779b97f4a7c16ULL + (seed << 6) + (seed >> 2);
				};

			hash_combine(h, v.position.x);
			hash_combine(h, v.position.y);
			hash_combine(h, v.position.z);
			hash_combine(h, v.uv.x);
			hash_combine(h, v.uv.y);
			hash_combine(h, v.normal.x);
			hash_combine(h, v.normal.y);
			hash_combine(h, v.normal.z);
			return h;
		}
	};

	struct Face {
		std::vector<uint32_t> vertexIndices{};
	};

	class Mesh {
	public:
		// Constructors/Destructors
		Mesh() = default;
		explicit Mesh(const std::string& filePath);
		~Mesh() = default;

		// Core functionality
		bool loadFromFile(const std::string& filePath);
		void clear();

		// Geometry access and manipulation
		const std::vector<Face>& getFaces() const { return faces; }
		const std::vector<Vertex>& getVertices() const { return vertices; }

		// File path management
		const std::string& getSourcePath() const { return sourcePath; }

	private:
		// Core data
		std::vector<Vertex> vertices;
		std::vector<Face> faces;
		std::string sourcePath;

		// Geometry processing
		bool parseGeometry(const std::string& fileContent);

		// Caching and optimization
		std::unordered_map<std::string, size_t> vertexCache;
		void buildVertexCache();
		void clearCache();
	};
}