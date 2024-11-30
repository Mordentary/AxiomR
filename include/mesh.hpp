#pragma once
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

#include"renderer.hpp"
#include"vec.hpp"

namespace AR {
	struct Vertex {
		Vec3f position;
		Vec3f normal;
		Vec3f uv;
	};

	struct Face {
		std::vector<uint32_t> vertexIndices;
		std::vector<uint32_t> normalIndices;
		std::vector<uint32_t> uvIndices;
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
		const std::vector<Vertex>& getVertices() const { return vertices; }
		const std::vector<Face>& getFaces() const { return faces; }

		// Geometry operations
		void computeNormals();
		void optimizeGeometry();
		bool isValid() const;

		// File path management
		const std::string& getSourcePath() const { return sourcePath; }

	private:
		// Core data
		std::vector<Vertex> vertices;
		std::vector<Face> faces;
		std::string sourcePath;

		// Geometry processing
		bool parseGeometry(const std::string& fileContent);
		void generateSmoothNormals();
		void generateFaceNormals();

		// Caching and optimization
		std::unordered_map<std::string, size_t> vertexCache;
		void buildVertexCache();
		void clearCache();
	};
}