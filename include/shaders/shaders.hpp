#include"../pipeline.hpp"
#include"IShader.hpp"
namespace AR
{
	struct FlatShader : public IShader {
		friend Pipeline;
	public:
		virtual ~FlatShader() {}

		virtual Vec4f vertex(const Vertex& vertex, int indexInsideFace) {
			Vec4f gl_Vertex = toVec4f(vertex.position, 1.0f);
			gl_Vertex = mvp * gl_Vertex;
			varying_normal.set_col(indexInsideFace, perspectiveDivision(gl_Vertex));
			return (gl_Vertex);
		}

		virtual bool fragment(Vec3f& bar, Vec4f& color)
		{
			Vec3f n = (varying_normal.get_col(1) - varying_normal.get_col(0)).cross(varying_normal.get_col(2) - varying_normal.get_col(0));
			n.normalize();
			float intensity = lightDirection.dot(n);
			intensity = std::clamp(intensity, 0.0f, 1.0f);
			color = Vec4f{ 1.0f, 1.0f, 1.0f, 1.0f } *intensity;;
			return false;
		}

	public:
		Vec3f lightDirection;
	private:
		Mat<3, 3, float> varying_normal;
	};

	struct DefaultShader : public IShader
	{
		friend Pipeline;
	public:
		virtual ~DefaultShader() {}

		virtual Vec4f vertex(const Vertex& vertex, int indexInsideFace) {
			Vec4f gl_Vertex = toVec4f(vertex.position, 1.0f);
			gl_Vertex = mvp * gl_Vertex;
			varying_uv.set_col(indexInsideFace, vertex.uv);
			Vec3f T = (toVec3f(model * toVec4f(vertex.tangent, 0.0)));
			T.normalize();
			Vec3f B = (toVec3f(model * toVec4f(vertex.bitangent, 0.0)));
			B.normalize();
			Vec3f N = (toVec3f(model * toVec4f(vertex.normal, 0.0)));
			N.normalize();

			TBN.set_col(0, T);
			TBN.set_col(1, B);
			TBN.set_col(2, N);
			TBN = TBN.transpose();
			Vec3f lightDirTangent = TBN * lightDirection;
			lightDirTangent.normalize();
			varying_lightDir.set_col(indexInsideFace, lightDirTangent);

			return (gl_Vertex);
		}

		virtual bool fragment(Vec3f& bar, Vec4f& color)
		{
			Vec2f uv = varying_uv * bar;
			Vec4f normMapValue = material->bumpTexture->sample(uv);
			Vec3f normal = toVec3f(normMapValue);
			normal = normal * 2.0f - 1.0f;
			normal.normalize();

			normal = -normal;
			Vec3f lightDir = varying_lightDir * bar;
			float intensity = lightDir.dot(normal);
			intensity = std::clamp(intensity, 0.0f, 1.0f);
			color = material->diffuseTexture->sample(uv) * intensity;
			//color = toVec4f(normal, 1.0f);
			return false;
		}
	public:
		Vec3f lightDirection;
	private:
		Mat<2, 3, float> varying_uv;
		Mat<3, 3, float> varying_lightDir;
		mat3f TBN;
	};

	struct PhongShader : public IShader {
		friend Pipeline;
	public:
		virtual ~PhongShader() {}

		virtual Vec4f vertex(const Vertex& vertex, int indexInsideFace) {
			Vec4f gl_Vertex = toVec4f(vertex.position, 1.0f);
			gl_Vertex = mvp * gl_Vertex;
			varying_uv.set_col(indexInsideFace, vertex.uv);
			Vec4f worldPos = model * toVec4f(vertex.position, 1.0f);
			varying_fragWorldPos.set_col(indexInsideFace, toVec3f(worldPos));

			mat3f worldTBN;
			Vec3f T = (toVec3f(model * toVec4f(vertex.tangent, 0.0))).normalized();
			Vec3f B = (toVec3f(model * toVec4f(vertex.bitangent, 0.0))).normalized();
			Vec3f N = (toVec3f(model * toVec4f(vertex.normal, 0.0))).normalized();
			worldTBN.set_col(0, T);
			worldTBN.set_col(1, B);
			worldTBN.set_col(2, N);

			// Pass the worldTBN to the appropriate varying based on the vertex index
			if (indexInsideFace == 0) varying_tbn0 = worldTBN;
			else if (indexInsideFace == 1) varying_tbn1 = worldTBN;
			else if (indexInsideFace == 2) varying_tbn2 = worldTBN;

			return (gl_Vertex);
		}

		virtual bool fragment(Vec3f& bar, Vec4f& color)
		{
			Vec2f uv = varying_uv * bar;
			Vec4f normMapValue = material->bumpTexture->sample(uv);
			Vec3f normalMapSample = toVec3f(normMapValue) * 2.0f - 1.0f;
			normalMapSample.normalize();

			// Interpolate the TBN matrices using barycentric coordinates
			mat3f interpolatedTBN = varying_tbn0 * bar.x + varying_tbn1 * bar.y + varying_tbn2 * bar.z;

			// Extract the interpolated basis vectors
			Vec3f interpolatedT = interpolatedTBN.get_col(0);
			Vec3f interpolatedB = interpolatedTBN.get_col(1);
			Vec3f interpolatedN = interpolatedTBN.get_col(2);

			// Normalize the interpolated normal
			interpolatedN.normalize();

			// Re-orthogonalize tangent and bitangent
			Vec3f interpolatedTNorm = (interpolatedT - interpolatedN * interpolatedN.dot(interpolatedT)).normalized();
			Vec3f interpolatedBNorm = interpolatedN.cross(interpolatedTNorm);

			// Reconstruct the orthonormal interpolated TBN matrix
			mat3f finalInterpolatedTBN;
			finalInterpolatedTBN.set_col(0, interpolatedTNorm);
			finalInterpolatedTBN.set_col(1, interpolatedBNorm);
			finalInterpolatedTBN.set_col(2, interpolatedN);

			// 1. Sample Textures
			Vec4f albedo = material->diffuseTexture->sample(uv);

			// 2. Normal Mapping
			// Transform the normal map sample from tangent to world space
			Vec3f normal = finalInterpolatedTBN * normalMapSample;
			normal.normalize();
			//normal = -normal;

			// 3. Lighting Setup
			Vec3f fragPos = varying_fragWorldPos * bar; // Interpolated world position
			Vec3f viewDir = (cameraPosition - fragPos).normalized();
			Vec3f lightDir = -lightDirection;

			// 4. Ambient Lighting
			float ambientStrength = 0.1f;
			Vec3f ambient = ambientStrength * (lightColor);

			// 6. Specular Lighting
			float specularStrength = 0.5f;
			Vec3f reflectDir = reflect(-lightDir, normal);
			float spec = std::pow(std::max(viewDir.dot(reflectDir), 0.0f), material->specularExponent * 50.0f);
			Vec3f specular = specularStrength * spec * (lightColor);

			// 5. Diffuse Lighting
			float diff = std::max(normal.dot(lightDir), 0.0f);
			Vec3f diffuse = diff * lightColor;
			Vec3f finalColor = (ambient + diffuse + specular) * toVec3f(albedo);

			color = toVec4f(finalColor, 1.0f);

			return false;
		}
	public:
		Vec3f lightDirection;
		Vec3f lightColor;
	private:
		Mat<2, 3, float> varying_uv;
		Mat<3, 3, float> varying_fragWorldPos;
		mat3f varying_tbn0;
		mat3f varying_tbn1;
		mat3f varying_tbn2;
	};
	struct PBRShader : public IShader {
		friend Pipeline;
	public:
		virtual ~PBRShader() {}

		virtual Vec4f vertex(const Vertex& vertex, int indexInsideFace) {
			Vec4f gl_Vertex = toVec4f(vertex.position, 1.0f);
			gl_Vertex = mvp * gl_Vertex;
			varying_uv.set_col(indexInsideFace, vertex.uv);
			Vec4f worldPos = model * toVec4f(vertex.position, 1.0f);
			varying_fragWorldPos.set_col(indexInsideFace, toVec3f(worldPos));

			// Calculate world-space TBN matrix for the current vertex
			mat3f worldTBN;
			Vec3f T = (toVec3f(model * toVec4f(vertex.tangent, 0.0))).normalized();
			Vec3f B = (toVec3f(model * toVec4f(vertex.bitangent, 0.0))).normalized();
			Vec3f N = (toVec3f(model * toVec4f(vertex.normal, 0.0))).normalized();
			worldTBN.set_col(0, T);
			worldTBN.set_col(1, B);
			worldTBN.set_col(2, N);

			// Pass the worldTBN to the appropriate varying based on the vertex index
			if (indexInsideFace == 0) varying_tbn0 = worldTBN;
			else if (indexInsideFace == 1) varying_tbn1 = worldTBN;
			else if (indexInsideFace == 2) varying_tbn2 = worldTBN;

			return (gl_Vertex);
		}

		virtual bool fragment(Vec3f& bar, Vec4f& color)
		{
			Vec2f uv = varying_uv * bar;
			Vec4f normMapValue = material->bumpTexture->sample(uv);
			Vec3f normalMapSample = toVec3f(normMapValue) * 2.0f - 1.0f;
			normalMapSample.normalize();

			// Interpolate the TBN matrices using barycentric coordinates
			mat3f interpolatedTBN = varying_tbn0 * bar.x + varying_tbn1 * bar.y + varying_tbn2 * bar.z;

			// Extract the interpolated basis vectors
			Vec3f interpolatedT = interpolatedTBN.get_col(0);
			Vec3f interpolatedB = interpolatedTBN.get_col(1);
			Vec3f interpolatedN = interpolatedTBN.get_col(2);

			// Normalize the interpolated normal
			interpolatedN.normalize();

			// Re-orthogonalize tangent and bitangent
			Vec3f interpolatedTNorm = (interpolatedT - interpolatedN * interpolatedN.dot(interpolatedT)).normalized();
			//Vec3f interpolatedBNorm = interpolatedN.cross(interpolatedTNorm);
			Vec3f interpolatedBNorm = (interpolatedB - interpolatedB.dot(interpolatedT) * interpolatedT - interpolatedB.dot(interpolatedN) * interpolatedN).normalized();
			interpolatedN = interpolatedTNorm.cross(interpolatedB).normalized();

			// Reconstruct the orthonormal interpolated TBN matrix
			mat3f finalInterpolatedTBN;
			finalInterpolatedTBN.set_col(0, interpolatedTNorm);
			finalInterpolatedTBN.set_col(1, interpolatedBNorm);
			finalInterpolatedTBN.set_col(2, interpolatedN);

			// 1. Sample Textures
			Vec3f albedo = toVec3f(material->diffuseTexture->sample(uv)).pow(Vec3f(2.2f));
			Vec4f metallicVec = material->metallicTexture->sample(uv);
			Vec4f roughnessVec = material->roughnessTexture ? material->roughnessTexture->sample(uv) : 1.0f;
			Vec4f aoVec = material->aoTexture ? material->aoTexture->sample(uv) : Vec4f(1.0f);

			float metallic = metallicVec.x;
			float roughness = roughnessVec.x;
			float ao = aoVec.x;
			// 2. Normal Mapping
			// Transform the normal map sample from tangent to world space
			Vec3f normal = finalInterpolatedTBN * normalMapSample;
			normal.normalize();
			//normal = -normal;

			// 3. Lighting Setup
			Vec3f fragPos = varying_fragWorldPos * bar; // Interpolated world position
			Vec3f viewDir = (cameraPosition - fragPos).normalized();
			Vec3f lightDir = -lightDirection;
			Vec3f halfwayDir = (lightDir + viewDir).normalized();

			// 4. PBR Calculations
			// a. Fresnel (Schlick's approximation)
			Vec3f F0 = lerp(Vec3f(0.04f), (albedo), metallic); // Base reflectivity (0.04 for dielectrics)
			Vec3f F = fresnelSchlick(std::max(halfwayDir.dot(normal), 0.0f), F0);

			// b. Normal Distribution Function (Trowbridge-Reitz GGX)
			float NDF = distributionGGX(normal, halfwayDir, roughness);

			// c. Geometry Function (Smith)
			float G = geometrySmith(normal, viewDir, lightDir, roughness);

			// d. Cook-Torrance BRDF
			Vec3f numerator = F * NDF * G;
			float denominator = 4.0f * std::max(normal.dot(viewDir), 0.0f) * std::max(normal.dot(lightDir), 0.0f) + 0.0001f; // Avoid division by zero
			Vec3f specular = numerator / denominator;

			// e. Diffuse (Lambertian)
			Vec3f kS = F;
			Vec3f kD = Vec3f(1.0f) - kS;
			kD *= (1.0f - metallic);

			// 5. Combine Lighting
			float NdotL = std::max(normal.dot(lightDir), 0.0f);
			Vec3f radiance = lightColor;
			Vec3f Lo = (kD * albedo / std::numbers::pi + specular) * radiance * NdotL;

			Vec3f ambient = Vec3f(0.03) * albedo * ao;
			Vec3f finalColor = ambient + Lo;
			//color = toVec4f((ambient + (diffuse + specular) * radiance), 1.0f);

			 //6. HDR and Tone Mapping (Reinhard)
			finalColor = finalColor / (finalColor + Vec3f(1.0f));
			//7. Gamma Correction
			color = toVec4f((finalColor).pow(Vec3f(1.0f / 2.2f)), 1.0f);

			return false;
		}
	public:
		Vec3f lightDirection;
		Vec3f lightColor;

	private:
		Vec3f fresnelSchlick(float cosTheta, Vec3f F0)
		{
			return F0 + (Vec3f(1.0f) - F0) * std::pow(std::clamp(1.0f - cosTheta, 0.0f, 1.0f), 5.0f);
		}

		// Normal Distribution Function (Trowbridge-Reitz GGX)
		float distributionGGX(Vec3f N, Vec3f H, float roughness)
		{
			float a = roughness * roughness;
			float a2 = a * a;
			float NdotH = std::max(N.dot(H), 0.0f);
			float NdotH2 = NdotH * NdotH;

			float num = a2;
			float denom = (NdotH2 * (a2 - 1.0f) + 1.0f);
			denom = std::numbers::pi * denom * denom;

			return num / denom;
		}

		// Geometry Function (Schlick-GGX)
		float geometrySchlickGGX(float NdotV, float roughness)
		{
			float r = (roughness + 1.0);
			float k = (r * r) / 8.0;

			float num = NdotV;
			float denom = NdotV * (1.0 - k) + k;
			return num / denom;
		}

		// Geometry Function (Smith)
		float geometrySmith(Vec3f N, Vec3f V, Vec3f L, float roughness)
		{
			float NdotV = std::max(N.dot(V), 0.0f);
			float NdotL = std::max(N.dot(L), 0.0f);
			float ggx1 = geometrySchlickGGX(NdotV, roughness);
			float ggx2 = geometrySchlickGGX(NdotL, roughness);

			return ggx1 * ggx2;
		}

		// Linear interpolation
		Vec3f lerp(Vec3f a, Vec3f b, float t)
		{
			return a + t * (b - a);
		}
	private:
		Mat<2, 3, float> varying_uv;
		Mat<3, 3, float> varying_fragWorldPos;
		mat3f varying_tbn0;
		mat3f varying_tbn1;
		mat3f varying_tbn2;
	};
}