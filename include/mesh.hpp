#pragma once
#include <vector>
#include <array>
#include <string>
#include <memory>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_set>
#include <functional>
#include"renderer.hpp"
#include"vec.hpp"

namespace AR {
	struct Vertex {
		Vec3f position;
		Vec2f uv;
		Vec3f normal;

		bool operator==(const Vertex& other) const {
			return position == other.position &&
				uv == other.uv &&
				normal == other.normal;
		}
	};

	struct Face {
		std::vector<uint32_t> vertexIndices{ 0,0,0 };
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

namespace std {
	template<>
	struct hash<AR::Vertex> {
		size_t operator()(const AR::Vertex& vertex) const {
			return ((hash<float>()(vertex.position.x) ^
				(hash<float>()(vertex.position.y) << 1)) >> 1) ^
				(hash<float>()(vertex.position.z) << 1) ^
				(hash<float>()(vertex.uv.x) << 2) ^
				(hash<float>()(vertex.uv.y) << 3) ^
				(hash<float>()(vertex.normal.x) << 4) ^
				(hash<float>()(vertex.normal.y) << 5) ^
				(hash<float>()(vertex.normal.z) << 6);
		}
	};
}