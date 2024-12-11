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

	bool Mesh::parseGeometry(const std::string& fileContent) {
		std::istringstream stream(fileContent);
		std::string line;

		std::vector<Vec3f> positions;
		std::vector<Vec2f> uvs;
		std::vector<Vec3f> normals;

		std::vector<Vertex> vertices;
		std::unordered_map<Vertex, uint32_t, VertexHash> uniqueVertices;

		std::vector<Face> faces;

		while (std::getline(stream, line)) {
			if (line.empty()) continue;
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
				uvs.emplace_back(u, v);
			}
			else if (type == "vn") {
				float nx, ny, nz;
				lineStream >> nx >> ny >> nz;
				normals.emplace_back(nx, ny, nz);
			}
			else if (type == "f") {
				Face face;
				std::string token;
				while (lineStream >> token) {
					int vIndex = -1, vtIndex = -1, vnIndex = -1;
					// Split token by '/'
					std::vector<std::string> parts;
					{
						std::stringstream ss(token);
						std::string sub;
						while (std::getline(ss, sub, '/'))
							parts.push_back(sub);
					}

					if (!parts.empty() && !parts[0].empty()) vIndex = std::stoi(parts[0]) - 1;
					if (parts.size() > 1 && !parts[1].empty()) vtIndex = std::stoi(parts[1]) - 1;
					if (parts.size() > 2 && !parts[2].empty()) vnIndex = std::stoi(parts[2]) - 1;

					if (vIndex < 0 || vIndex >= (int)positions.size()) continue;

					Vertex vertexData;
					vertexData.position = positions[vIndex];
					vertexData.uv = (vtIndex >= 0 && vtIndex < (int)uvs.size()) ? uvs[vtIndex] : Vec2f(0.0f, 0.0f);
					vertexData.normal = (vnIndex >= 0 && vnIndex < (int)normals.size()) ? normals[vnIndex] : Vec3f(0.0f, 0.0f, 0.0f);

					uint32_t vertexIndex;
					if (uniqueVertices.find(vertexData) == uniqueVertices.end()) {
						vertexIndex = (uint32_t)vertices.size();
						uniqueVertices[vertexData] = vertexIndex;
						vertices.push_back(vertexData);
					}
					else {
						vertexIndex = uniqueVertices[vertexData];
					}

					face.vertexIndices.push_back(vertexIndex);
				}
				faces.push_back(face);
			}
		}

		this->vertices = vertices;
		this->faces = faces;

		return true;
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
}