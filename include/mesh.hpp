#pragma once
#include <vector>
#include <unordered_map>
#include <string>
#include"math.hpp"
#include"texture.hpp"

namespace AR {
	struct Vertex {
		glm::vec3 position;
		glm::vec2 uv;
		glm::vec3 normal;
		glm::vec3 tangent;
		glm::vec3 bitangent;
		bool operator==(const Vertex& other) const {
			return position == other.position && uv == other.uv && normal == other.normal;
		}
	};

	struct Material {
		std::string name;
		std::unique_ptr<Texture> diffuseTexture;
		std::unique_ptr<Texture> metallicTexture;
		std::unique_ptr<Texture> bumpTexture;
		std::unique_ptr<Texture> roughnessTexture;
		std::unique_ptr<Texture> aoTexture;
		glm::vec3 ambient;
		glm::vec3 diffuse;
		glm::vec3 specular;
		float specularExponent;
		float refractiveIndex;
		float dissolve;
		int illuminationModel;
	};
	struct MaterialGroup {
		std::string materialName;
		size_t startIndex;
		size_t faceCount;
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
		Mesh() = default;
		explicit Mesh(const std::string& filePath);
		~Mesh() = default;

		bool loadFromFile(const std::string& filePath);
		void clear();

		const std::vector<Face>& getFaces() const { return m_Faces; }
		const std::vector<Vertex>& getVertices() const { return m_Vertices; }
		const std::string& getSourcePath() const { return m_ModelSrcPath; }

		const Material* getMaterial(const std::string& matName) const {
			return m_Materials.at(matName).get();
		}

		const std::vector<MaterialGroup>& getMaterialGroups() const {
			return m_MaterialGroups;
		}

	private:
		std::vector<Vertex> m_Vertices;
		std::vector<Face> m_Faces;
		std::string m_ModelSrcPath;
		std::unordered_map<std::string, std::unique_ptr<Material>> m_Materials;
		std::vector<MaterialGroup> m_MaterialGroups;

		void calculateTangentBitangent();
		bool parseModelFile(const std::string& fileContent);
		bool loadMaterial();
		void parseMaterialData(const std::string& type, std::istringstream& lineStream, Material& currentMaterial, const std::string& mtlFilePath);
		std::string getBaseName(const std::string& filePath);
		std::string getDirectory(const std::string& filePath);
	};
}