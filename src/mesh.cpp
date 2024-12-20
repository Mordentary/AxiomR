#include "mesh.hpp"
#include "texture.hpp"
#include <fstream>
#include <sstream>
#include <map>

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
		m_ModelSrcPath = filePath;
		return parseModelFile(buffer.str());
	}

	void Mesh::clear() {
		m_Vertices.clear();
		m_Faces.clear();
		m_Materials.clear();
	}

	std::string Mesh::getBaseName(const std::string& filePath) {
		size_t lastSlash = filePath.find_last_of("/\\");
		size_t lastDot = filePath.find_last_of(".");

		if (lastDot == std::string::npos) {
			return (lastSlash == std::string::npos) ? filePath : filePath.substr(lastSlash + 1);
		}

		size_t start = (lastSlash == std::string::npos) ? 0 : lastSlash + 1;
		return filePath.substr(start, lastDot - start);
	}

	// Helper function to extract the directory from a file path
	std::string Mesh::getDirectory(const std::string& filePath) {
		size_t lastSlash = filePath.find_last_of("/\\");
		if (lastSlash != std::string::npos) {
			return filePath.substr(0, lastSlash);
		}
		return ""; // Current directory if no slash found
	}

	void processTexturePathLine(std::istringstream& lineStream, std::string& texturePath)
	{
		std::getline(lineStream, texturePath);
		size_t first = texturePath.find_first_not_of(" \t");
		if (first == std::string::npos) {
			std::cerr << "Error: Empty texture path." << std::endl;
		}
		size_t last = texturePath.find_last_not_of(" \t");
		texturePath = texturePath.substr(first, (last - first + 1));
	}

	void Mesh::parseMaterialData(const std::string& type, std::istringstream& lineStream, Material& currentMaterial, const std::string& mtlFilePath) {
		if (type == "Ka")
		{
			float r, g, b;
			lineStream >> r >> g >> b;
			currentMaterial.ambient = Vec3f{ r,g,b };
		}
		else if (type == "Kd")
		{
			float r, g, b;
			lineStream >> r >> g >> b;
			currentMaterial.diffuse = Vec3f{ r,g,b };
		}
		else if (type == "Ks")
		{
			float r, g, b;
			lineStream >> r >> g >> b;
			currentMaterial.specular = Vec3f{ r,g,b };
		}
		else if (type == "Ns")
		{
			lineStream >> currentMaterial.specularExponent;
		}
		else if (type == "Ni")
		{
			lineStream >> currentMaterial.refractiveIndex;
		}
		else if (type == "d")
		{
			lineStream >> currentMaterial.dissolve;
		}
		else if (type == "illum")
		{
			lineStream >> currentMaterial.illuminationModel;
		}
		else if (type == "map_Kd")
		{
			std::string texturePath;
			processTexturePathLine(lineStream, texturePath);
			std::string fullTexturePath = getDirectory(mtlFilePath) + "/" + texturePath;
			currentMaterial.diffuseTexture = std::make_unique<Texture>(fullTexturePath);
		}
		else if (type == "map_Ks")
		{
			std::string texturePath;
			processTexturePathLine(lineStream, texturePath);
			std::string fullTexturePath = getDirectory(mtlFilePath) + "/" + texturePath;
			currentMaterial.metallicTexture = std::make_unique<Texture>(fullTexturePath);
		}
		else if (type == "map_Ns")
		{
			std::string texturePath;
			processTexturePathLine(lineStream, texturePath);
			std::string fullTexturePath = getDirectory(mtlFilePath) + "/" + texturePath;
			currentMaterial.roughnessTexture = std::make_unique<Texture>(fullTexturePath);
		}
		else if (type == "map_Bump" || type == "bump" || type == "norm")
		{
			std::string texturePath;
			processTexturePathLine(lineStream, texturePath);
			std::string fullTexturePath = getDirectory(mtlFilePath) + "/" + texturePath;
			currentMaterial.bumpTexture = std::make_unique<Texture>(fullTexturePath);
		}
		else if (type == "map_A0")
		{
			std::string texturePath;
			processTexturePathLine(lineStream, texturePath);
			std::string fullTexturePath = getDirectory(mtlFilePath) + "/" + texturePath;
			currentMaterial.aoTexture = std::make_unique<Texture>(fullTexturePath);
		}
		//else if (type == "refl")
		//{
		//	std::string texturePath;
		//	lineStream >> texturePath;
		//	std::string fullTexturePath = getDirectory(mtlFilePath) + "/" + texturePath;
		//	currentMaterial.reflectionTexture = std::make_unique<Texture>(fullTexturePath);
		//}
	}
	uint32_t countMaterialsInMtlFile(const std::string& filePath) {
		std::ifstream file(filePath);
		if (!file.is_open()) {
			return -1; // Indicate error if file can't be opened
		}

		std::string line;
		int materialCount = 0;
		bool foundCount = false;

		while (std::getline(file, line)) {
			std::istringstream lineStream(line);
			std::string token;
			lineStream >> token;
			if (token.empty()) continue; //skip blank lines

			if (token == "#") {
				lineStream >> token;
				if (token == "Material") {
					lineStream >> token;
					if (token == "Count:")
					{
						lineStream >> materialCount;
						foundCount = true;
					}
				}
			}
			else if (token == "newmtl") {
				if (!foundCount)
					materialCount++;
			}
		}

		file.close();
		return materialCount;
	}

	bool Mesh::loadMaterial() {
		std::string mtlFilePath = getDirectory(m_ModelSrcPath) + "/" + getBaseName(m_ModelSrcPath) + ".mtl";
		uint32_t matCount = countMaterialsInMtlFile(mtlFilePath);
		int matIndex = -1;
		std::ifstream file(mtlFilePath);
		if (!file.is_open()) {
			return false;
		}

		std::string line;
		m_Materials.reserve(matCount);
		std::string name;
		while (std::getline(file, line)) {
			if (line.empty()) continue;

			std::istringstream lineStream(line);
			std::string type;
			lineStream >> type;

			if (type == "newmtl")
			{
				if (matIndex < int(matCount - 1))
					++matIndex;

				name.clear();
				lineStream >> name;
				m_Materials[name] = std::make_unique<Material>();
				lineStream >> m_Materials[name]->name;
			}
			else
			{
				if (matIndex != -1 && m_Materials[name])

					parseMaterialData(type, lineStream, *m_Materials[name], mtlFilePath);
			}
		}

		file.close();
		return true;
	}
	void Mesh::calculateTangentBitangent() {
		std::map<uint32_t, Vec3f> tangentMap;
		std::map<uint32_t, Vec3f> bitangentMap;

		for (const auto& face : m_Faces) {
			if (face.vertexIndices.size() < 3) continue;

			const Vertex& v0 = m_Vertices[face.vertexIndices[0]];
			const Vertex& v1 = m_Vertices[face.vertexIndices[1]];
			const Vertex& v2 = m_Vertices[face.vertexIndices[2]];

			Vec3f edge1 = v1.position - v0.position;
			Vec3f edge2 = v2.position - v0.position;

			Vec2f deltaUV1 = v1.uv - v0.uv;
			Vec2f deltaUV2 = v2.uv - v0.uv;

			float denominator = (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);
			if (fabs(denominator) < 1e-8f) {
				Vec3f fallbackTangent = (fabs(v0.normal.x) > 0.8f) ? Vec3f(0.0f, 1.0f, 0.0f) : Vec3f(1.0f, 0.0f, 0.0f);
				Vec3f fallbackBitangent = v0.normal.cross(fallbackTangent);

				for (uint32_t index : face.vertexIndices) {
					tangentMap[index] += fallbackTangent;
					bitangentMap[index] += fallbackBitangent;
				}
				continue;
			}

			float f = 1.0f / denominator;
			Vec3f tangent = f * (deltaUV2.y * edge1 - deltaUV1.y * edge2);
			Vec3f bitangent = f * (-deltaUV2.x * edge1 + deltaUV1.x * edge2);

			Vec3f expectedBitangent = v0.normal.cross(tangent);
			if (expectedBitangent.dot(bitangent) < 0.0f) {
				bitangent = -bitangent;
			}

			for (uint32_t index : face.vertexIndices) {
				tangentMap[index] += tangent;
				bitangentMap[index] += bitangent;
			}
		}

		for (auto& vertex : m_Vertices) {
			uint32_t index = static_cast<uint32_t>(&vertex - &m_Vertices[0]);

			Vec3f normal = vertex.normal;
			Vec3f tangent = (tangentMap.count(index) > 0) ? tangentMap[index] : Vec3f(1.0f, 0.0f, 0.0f);
			Vec3f bitangent = (bitangentMap.count(index) > 0) ? bitangentMap[index] : Vec3f(0.0f, 1.0f, 0.0f);

			// 1. Orthonormalize: Make tangent perpendicular to normal
			tangent = (tangent - normal * normal.dot(tangent));
			tangent.normalize();
			// 2. Recompute bitangent using cross product: This ensures it's perpendicular to both normal and tangent
			bitangent = normal.cross(tangent);

			// 3. Handedness check:
			float handedness = (bitangent.dot(bitangentMap[index]) < 0.0f) ? -1.0f : 1.0f;
			bitangent *= handedness;

			if (tangent.length() < 1e-8f) {
				tangent = (fabs(normal.x) > 0.8f) ? Vec3f(0.0f, 1.0f, 0.0f) : Vec3f(1.0f, 0.0f, 0.0f);
				tangent = ((tangent - normal) * normal.dot(tangent));
				tangent.normalize();
			}

			if (bitangent.length() < 1e-8f) {
				bitangent = (fabs(normal.z) > 0.8f) ? Vec3f(1.0f, 0.0f, 0.0f) : Vec3f(0.0f, 0.0f, 1.0f);
				bitangent = normal.cross(tangent);
			}

			// Assign final values
			vertex.tangent = tangent;
			vertex.bitangent = bitangent;
		}
	}

	bool Mesh::parseModelFile(const std::string& fileContent) {
		std::istringstream stream(fileContent);
		std::string line;

		std::vector<Vec3f> positions;
		std::vector<Vec2f> uvs;
		std::vector<Vec3f> normals;

		std::vector<Vertex> vertices;
		std::unordered_map<Vertex, uint32_t, VertexHash> uniqueVertices;

		std::vector<Face> faces;
		std::string currentMaterialName = "";
		size_t currentGroupStart = 0; // Start index of the current group
		m_MaterialGroups.clear();    // Clear previous material groups

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
			else if (type == "usemtl") {
				// Material change
				lineStream >> currentMaterialName;

				// If not the first group, add the previous group to the list
				if (!m_MaterialGroups.empty()) {
					m_MaterialGroups.back().faceCount = faces.size() - currentGroupStart;
				}

				// Add a new group with the current material (faceCount will be updated later)
				m_MaterialGroups.push_back({ currentMaterialName, faces.size(), 0 });
				currentGroupStart = faces.size();
			}
			else if (type == "f") {
				std::vector<uint32_t> vertexIndices; // Store vertex indices for the current face
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

					vertexIndices.push_back(vertexIndex); // Add index to the current face's list
				}

				// Triangulate the face (if it has more than 3 vertices)
				if (vertexIndices.size() >= 3) {
					for (size_t i = 1; i + 1 < vertexIndices.size(); ++i) {
						Face newFace;
						newFace.vertexIndices.push_back(vertexIndices[0]);     // First vertex
						newFace.vertexIndices.push_back(vertexIndices[i]);     // Second vertex
						newFace.vertexIndices.push_back(vertexIndices[i + 1]); // Third vertex
						faces.push_back(newFace);
					}
				}
			}
			if (!m_MaterialGroups.empty()) {
				m_MaterialGroups.back().faceCount = faces.size() - currentGroupStart;
			}
		}

		this->m_Vertices = vertices;
		this->m_Faces = faces;

		if (!loadMaterial())
		{
			std::cerr << "Warning: Could not load material file: " << getDirectory(m_ModelSrcPath) + "/" + getBaseName(m_ModelSrcPath) + ".mtl" << std::endl;
		}
		calculateTangentBitangent();

		return true;
	}
}