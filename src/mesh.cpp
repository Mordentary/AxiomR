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
		std::unordered_map<Vertex, uint32_t> uniqueVertices;

		std::vector<Face> faces; // Store faces here

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
				uvs.emplace_back(u, v);
			}
			else if (type == "vn") {
				float nx, ny, nz;
				lineStream >> nx >> ny >> nz;
				normals.emplace_back(nx, ny, nz);
			}
			else if (type == "f") {
				Face face;

				std::string vertexStr;
				while (lineStream >> vertexStr) {
					int vIndex = -1, vtIndex = -1, vnIndex = -1;

					size_t firstSlash = vertexStr.find('/');
					size_t secondSlash = vertexStr.find('/', firstSlash + 1);

					try {
						// Only vertex index
						if (firstSlash == std::string::npos) {
							vIndex = std::stoi(vertexStr) - 1;
						}
						else {
							// Extract vertex index
							std::string vIndexStr = vertexStr.substr(0, firstSlash);
							vIndex = !vIndexStr.empty() ? std::stoi(vIndexStr) - 1 : -1;

							// Check for double slash (v//vn)
							if (vertexStr[firstSlash + 1] == '/') {
								// v//vn format
								size_t vnStart = firstSlash + 2;
								std::string vnIndexStr = vertexStr.substr(vnStart);
								vnIndex = !vnIndexStr.empty() ? std::stoi(vnIndexStr) - 1 : -1;
							}
							else {
								// v/vt or v/vt/vn format
								size_t vtStart = firstSlash + 1;
								size_t vtEnd = (secondSlash == std::string::npos) ? std::string::npos : secondSlash - vtStart;
								std::string vtIndexStr = vertexStr.substr(vtStart, vtEnd);
								vtIndex = !vtIndexStr.empty() ? std::stoi(vtIndexStr) - 1 : -1;

								if (secondSlash != std::string::npos) {
									std::string vnIndexStr = vertexStr.substr(secondSlash + 1);
									vnIndex = !vnIndexStr.empty() ? std::stoi(vnIndexStr) - 1 : -1;
								}
							}
						}
					}
					catch (const std::exception& e) {
						// Handle error: invalid index
						continue; // Skip this vertex
					}

					// Now, get the vertex data
					Vertex vertexData = {};

					// Get position
					if (vIndex >= 0 && vIndex < positions.size()) {
						vertexData.position = positions[vIndex];
					}
					else {
						continue; // Skip this vertex
					}

					if (vtIndex >= 0 && vtIndex < uvs.size()) {
						vertexData.uv = uvs[vtIndex];
					}
					else {
						vertexData.uv = Vec2f(0.0f, 0.0f); // Default UV
					}

					// Get normal
					if (vnIndex >= 0 && vnIndex < normals.size()) {
						vertexData.normal = normals[vnIndex];
					}
					else {
						vertexData.normal = Vec3f(0.0f, 0.0f, 0.0f); // Default normal
					}

					uint32_t vertexIndex;
					if (uniqueVertices.count(vertexData) == 0) {
						vertexIndex = static_cast<uint32_t>(vertices.size());
						uniqueVertices[vertexData] = vertexIndex;
						vertices.push_back(vertexData);
					}
					else {
						vertexIndex = uniqueVertices[vertexData];
					}

					// Corrected: Use push_back instead of indexing
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