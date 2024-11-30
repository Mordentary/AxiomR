#include "Mesh.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>

namespace AR {
	Mesh::Mesh(const std::string& filePath) {
		if (!loadFromFile(filePath)) {
			throw std::runtime_error("Failed to load mesh from file: " + filePath);
		}
	}

	bool Mesh::loadFromFile(const std::string& filePath) {
		std::ifstream file(filePath, std::ios::binary);
		if (!file.is_open()) {
			return false;
		}

		std::stringstream buffer;
		buffer << file.rdbuf();
		file.close();

		clear();
		sourcePath = filePath;
		return parseGeometry(buffer.str());
	}

	void Mesh::clear() {
		vertices.clear();
		faces.clear();
		clearCache();
	}

	void Mesh::computeNormals() {
		if (faces.size() > 0 && faces[0].normalIndices.empty()) {
			generateFaceNormals();
		}
		else {
			generateSmoothNormals();
		}
	}

	void Mesh::optimizeGeometry() {
		if (vertices.empty() || faces.empty()) {
			return;
		}

		buildVertexCache();

		std::vector<Vertex> optimizedVertices;
		std::vector<Face> optimizedFaces = faces;
		std::vector<uint32_t> vertexRemap(vertices.size());

		optimizedVertices.reserve(vertices.size());
		size_t currentIndex = 0;

		for (auto& face : optimizedFaces) {
			for (auto& index : face.vertexIndices) {
				if (vertexRemap[index] == 0 && index != 0) {
					optimizedVertices.push_back(vertices[index]);
					vertexRemap[index] = currentIndex++;
				}
				index = vertexRemap[index];
			}
		}

		vertices = std::move(optimizedVertices);
		faces = std::move(optimizedFaces);
	}

	bool Mesh::isValid() const {
		if (vertices.empty() || faces.empty()) {
			return false;
		}

		for (const auto& face : faces) {
			if (face.vertexIndices.size() < 3) {
				return false;
			}

			for (uint32_t index : face.vertexIndices) {
				if (index >= vertices.size()) {
					return false;
				}
			}
		}

		return true;
	}

	bool Mesh::parseGeometry(const std::string& fileContent) {
		std::istringstream stream(fileContent);
		std::string line;
		std::vector<Vec3f> positions;
		std::vector<Vec3f> uvs;
		std::vector<Vec3f> normals;

		while (std::getline(stream, line)) {
			std::istringstream lineStream(line);
			std::string type;
			lineStream >> type;

			if (type == "v") {
				float x, y, z;
				lineStream >> x >> y >> z;
				positions.emplace_back(x, y, z);
			}
			else if (type == "vt") {
				float u, v;
				lineStream >> u >> v;
				uvs.emplace_back(u, v, 0.0f);
			}
			else if (type == "vn") {
				float nx, ny, nz;
				lineStream >> nx >> ny >> nz;
				normals.emplace_back(nx, ny, nz);
			}
			else if (type == "f") {
				Face face;
				std::string vertex;
				while (lineStream >> vertex) {
					std::istringstream vertexStream(vertex);
					std::string indexStr;
					std::vector<uint32_t> indices;

					while (std::getline(vertexStream, indexStr, '/')) {
						if (!indexStr.empty()) {
							indices.push_back(std::stoul(indexStr) - 1);
						}
					}

					if (!indices.empty()) {
						face.vertexIndices.push_back(indices[0]);
						if (indices.size() > 1) {
							face.uvIndices.push_back(indices[1]);
						}
						if (indices.size() > 2) {
							face.normalIndices.push_back(indices[2]);
						}
					}
				}
				faces.push_back(face);
			}
		}

		for (size_t i = 0; i < positions.size(); ++i) {
			Vertex vertex;
			vertex.position = positions[i];

			if (i < normals.size()) {
				vertex.normal = normals[i];
			}

			if (i < uvs.size()) {
				vertex.uv = uvs[i];
			}

			vertices.push_back(vertex);
		}

		return isValid();
	}

	void Mesh::generateSmoothNormals() {
		std::vector<Vec3f> vertexNormals(vertices.size(), Vec3f(0, 0, 0));

		for (const auto& face : faces) {
			if (face.vertexIndices.size() < 3) continue;

			const Vec3f& v0 = vertices[face.vertexIndices[0]].position;
			const Vec3f& v1 = vertices[face.vertexIndices[1]].position;
			const Vec3f& v2 = vertices[face.vertexIndices[2]].position;

			Vec3f normal = (v1 - v0).cross(v2 - v0).normalized();

			for (uint32_t index : face.vertexIndices) {
				vertexNormals[index] = vertexNormals[index] + normal;
			}
		}

		for (size_t i = 0; i < vertices.size(); ++i) {
			if (vertexNormals[i].lengthSquared() > 0) {
				vertices[i].normal = vertexNormals[i].normalized();
			}
		}
	}

	void Mesh::generateFaceNormals() {
		for (auto& face : faces) {
			if (face.vertexIndices.size() < 3) continue;

			const Vec3f& v0 = vertices[face.vertexIndices[0]].position;
			const Vec3f& v1 = vertices[face.vertexIndices[1]].position;
			const Vec3f& v2 = vertices[face.vertexIndices[2]].position;

			Vec3f normal = (v1 - v0).cross(v2 - v0).normalized();

			for (uint32_t index : face.vertexIndices) {
				vertices[index].normal = normal;
			}
		}
	}

	void Mesh::buildVertexCache() {
		vertexCache.clear();
		for (size_t i = 0; i < vertices.size(); ++i) {
			const auto& v = vertices[i];
			std::string key = std::to_string(v.position.x) + "," +
				std::to_string(v.position.y) + "," +
				std::to_string(v.position.z);
			vertexCache[key] = i;
		}
	}

	void Mesh::clearCache() {
		vertexCache.clear();
	}
} // namespace AR