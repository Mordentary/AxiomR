#pragma once
#include "../pipeline.hpp"
#include "IShader.hpp"
#include "tracy/Tracy.hpp"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/component_wise.hpp>    // For glm::pow(vec, vec)
#include <glm/gtx/norm.hpp>             // For glm::length2(), etc.
#include <algorithm>                    // For std::clamp
#include <cmath>



namespace AR
{
	// --------------------------------------------------------------------------
	// FLAT SHADER
	// --------------------------------------------------------------------------
	struct FlatShader : public IShader
	{
		friend Pipeline;

	public:
		virtual ~FlatShader() {}

		virtual glm::vec4 vertex(const Vertex& vertex, int indexInsideFace) override
		{
			// Model-View-Projection transform
			glm::vec4 gl_Vertex = mvp * glm::vec4(vertex.position, 1.0f);

			// Store the normal in "varying_normal" in object/world space
			// for interpolation:
			//
			// We apply inverse-transpose(model) to the normal:
			glm::vec3 transformedNormal =
				glm::mat3(glm::transpose(glm::inverse(model))) * vertex.normal;
			varying_normal[indexInsideFace] = transformedNormal;

			return gl_Vertex;
		}

		virtual bool fragment(glm::vec3& bar, glm::vec4& color) override
		{
			// Interpolate the normal using the barycentric coordinates
			glm::vec3 n = bar.x * varying_normal[0]
				+ bar.y * varying_normal[1]
				+ bar.z * varying_normal[2];
			n = glm::normalize(n);

			// Basic lambertian-like intensity
			float intensity = std::clamp(glm::dot(-lightDirection, n), 0.0f, 1.0f);

			color = glm::vec4(1.0f) * intensity;
			return false;
		}

	public:
		glm::vec3 lightDirection;

	private:
		// For each face, we store the 3 normals in an array:
		glm::vec3 varying_normal[3];
	};

	// --------------------------------------------------------------------------
	// DEFAULT SHADER
	// --------------------------------------------------------------------------
	struct DefaultShader : public IShader
	{
		friend Pipeline;

	public:
		virtual ~DefaultShader() {}

		virtual glm::vec4 vertex(const Vertex& vertex, int indexInsideFace) override
		{
			// Transform position
			glm::vec4 gl_Vertex = mvp * glm::vec4(vertex.position, 1.0f);

			// Store UV
			varying_uv[indexInsideFace] = vertex.uv;

			// Transform tangent, bitangent, normal
			glm::vec3 T = glm::normalize(glm::vec3(model * glm::vec4(vertex.tangent, 0.0f)));
			glm::vec3 B = glm::normalize(glm::vec3(model * glm::vec4(vertex.bitangent, 0.0f)));
			glm::vec3 N = glm::normalize(glm::vec3(model * glm::vec4(vertex.normal, 0.0f)));

			// Build TBN; Note we transpose afterwards in original code
			glm::mat3 TBN;
			TBN[0] = T; // first column
			TBN[1] = B; // second column
			TBN[2] = N; // third column

			TBN = glm::transpose(TBN);

			// Light direction in tangent space
			glm::vec3 lightDirTangent = glm::normalize(TBN * lightDirection);
			varying_lightDir[indexInsideFace] = lightDirTangent;

			return gl_Vertex;
		}

		virtual bool fragment(glm::vec3& bar, glm::vec4& color) override
		{
			// Interpolate UV
			glm::vec2 uv = bar.x * varying_uv[0] +
				bar.y * varying_uv[1] +
				bar.z * varying_uv[2];

			// Sample normal map (in tangent space)
			glm::vec4 normMapValue = material->bumpTexture->sample(uv);
			glm::vec3 normal = glm::vec3(normMapValue) * 2.0f - glm::vec3(1.0f);
			normal = glm::normalize(normal);

			// Interpolate light direction (in tangent space)
			glm::vec3 lightDir = -(bar.x * varying_lightDir[0] +
				bar.y * varying_lightDir[1] +
				bar.z * varying_lightDir[2]);

			// Simple lambertian-like shading
			float intensity = glm::dot(lightDir, normal);
			intensity = std::clamp(intensity, 0.0f, 1.0f);

			// Sample diffuse
			glm::vec4 diffuseColor = material->diffuseTexture->sample(uv);
			color = diffuseColor * intensity;
			return false;
		}

	public:
		glm::vec3 lightDirection;

	private:
		// For each face, store 3 UV coords and 3 light directions
		glm::vec2 varying_uv[3];
		glm::vec3 varying_lightDir[3];
	};

	// --------------------------------------------------------------------------
	// PHONG SHADER
	// --------------------------------------------------------------------------
	struct PhongShader : public IShader
	{
		friend Pipeline;

	public:
		virtual ~PhongShader() {}

		virtual glm::vec4 vertex(const Vertex& vertex, int indexInsideFace) override
		{
			// Position -> clip space
			glm::vec4 gl_Vertex = mvp * glm::vec4(vertex.position, 1.0f);

			// Store UV
			varying_uv[indexInsideFace] = vertex.uv;

			// Compute and store world pos
			glm::vec4 worldPos = model * glm::vec4(vertex.position, 1.0f);
			varying_fragWorldPos[indexInsideFace] = glm::vec3(worldPos);

			// Build TBN in world space
			glm::vec3 T = glm::normalize(glm::vec3(model * glm::vec4(vertex.tangent, 0.0f)));
			glm::vec3 B = glm::normalize(glm::vec3(model * glm::vec4(vertex.bitangent, 0.0f)));
			glm::vec3 N = glm::normalize(glm::vec3(model * glm::vec4(vertex.normal, 0.0f)));

			glm::mat3 worldTBN;
			worldTBN[0] = T;
			worldTBN[1] = B;
			worldTBN[2] = N;

			// Store TBN in the appropriate slot
			if (indexInsideFace == 0) varying_tbn0 = worldTBN;
			else if (indexInsideFace == 1) varying_tbn1 = worldTBN;
			else if (indexInsideFace == 2) varying_tbn2 = worldTBN;

			return gl_Vertex;
		}

		virtual bool fragment(glm::vec3& bar, glm::vec4& color) override
		{
			// Interpolate UV
			glm::vec2 uv = bar.x * varying_uv[0] +
				bar.y * varying_uv[1] +
				bar.z * varying_uv[2];

			// Sample normal map (range [0,1] -> [-1,1])
			glm::vec4 normMapValue = material->bumpTexture->sample(uv);
			glm::vec3 normalMapSample = glm::vec3(normMapValue) * 2.0f - glm::vec3(1.0f);
			normalMapSample = glm::normalize(normalMapSample);

			// Interpolate TBN
			glm::mat3 interpolatedTBN =
				bar.x * varying_tbn0
				+ bar.y * varying_tbn1
				+ bar.z * varying_tbn2;

			// Extract T, B, N
			glm::vec3 interpolatedT = interpolatedTBN[0];
			glm::vec3 interpolatedB = interpolatedTBN[1];
			glm::vec3 interpolatedN = interpolatedTBN[2];

			// Normalize
			interpolatedN = glm::normalize(interpolatedN);

			// Re-orthogonalize T, B
			// T' = T - (N dot T)*N
			glm::vec3 interpolatedTNorm = glm::normalize(
				interpolatedT - interpolatedN * glm::dot(interpolatedN, interpolatedT)
			);
			glm::vec3 interpolatedBNorm = glm::cross(interpolatedN, interpolatedTNorm);

			// Reconstruct TBN
			glm::mat3 finalInterpolatedTBN;
			finalInterpolatedTBN[0] = interpolatedTNorm;
			finalInterpolatedTBN[1] = interpolatedBNorm;
			finalInterpolatedTBN[2] = interpolatedN;

			// Albedo from texture
			glm::vec4 albedo = material->diffuseTexture->sample(uv);

			// Transform normal from tangent to world
			glm::vec3 normal = glm::normalize(finalInterpolatedTBN * normalMapSample);

			// Lighting
			glm::vec3 fragPos = bar.x * varying_fragWorldPos[0]
				+ bar.y * varying_fragWorldPos[1]
				+ bar.z * varying_fragWorldPos[2];
			glm::vec3 viewDir = glm::normalize(cameraPosition - fragPos);
			glm::vec3 lightDir = -lightDirection;

			// Ambient
			float ambientStrength = 0.1f;
			glm::vec3 ambient = ambientStrength * lightColor;

			// Diffuse
			float diff = std::max(glm::dot(normal, lightDir), 0.0f);
			glm::vec3 diffuse = diff * lightColor;

			// Specular
			float specularStrength = 0.5f;
			glm::vec3 reflectDir = glm::reflect(-lightDir, normal);
			float spec = std::pow(std::max(glm::dot(viewDir, reflectDir), 0.0f),
				material->specularExponent * 50.0f);
			glm::vec3 specular = specularStrength * spec * lightColor;

			glm::vec3 finalColor = (ambient + diffuse + specular) * glm::vec3(albedo);

			color = glm::vec4(finalColor, 1.0f);
			return false;
		}

	public:
		glm::vec3 lightDirection;
		glm::vec3 lightColor;

	private:
		// Arrays to store each face's data
		glm::vec2 varying_uv[3];
		glm::vec3 varying_fragWorldPos[3];

		glm::mat3 varying_tbn0;
		glm::mat3 varying_tbn1;
		glm::mat3 varying_tbn2;
	};

	// --------------------------------------------------------------------------
	// PBR SHADER
	// --------------------------------------------------------------------------
	struct PBRShader : public IShader
	{
		friend Pipeline;

	public:
		virtual ~PBRShader() {}

		virtual glm::vec4 vertex(const Vertex& vertex, int indexInsideFace) override
		{
			ZoneScoped;

			glm::vec4 gl_Vertex = mvp * glm::vec4(vertex.position, 1.0f);

			// Store UV
			varying_uv[indexInsideFace] = vertex.uv;

			// Store world position
			glm::vec4 worldPos = model * glm::vec4(vertex.position, 1.0f);
			varying_fragWorldPos[indexInsideFace] = glm::vec3(worldPos);

			// Build TBN in world space
			glm::vec3 T = glm::normalize(glm::vec3(model * glm::vec4(vertex.tangent, 0.0f)));
			glm::vec3 B = glm::normalize(glm::vec3(model * glm::vec4(vertex.bitangent, 0.0f)));
			glm::vec3 N = glm::normalize(glm::vec3(model * glm::vec4(vertex.normal, 0.0f)));

			glm::mat3 worldTBN;
			worldTBN[0] = T;
			worldTBN[1] = B;
			worldTBN[2] = N;

			// Pass TBN to appropriate varying
			if (indexInsideFace == 0) varying_tbn0 = worldTBN;
			else if (indexInsideFace == 1) varying_tbn1 = worldTBN;
			else if (indexInsideFace == 2) varying_tbn2 = worldTBN;

			return gl_Vertex;
		}

		virtual bool fragment(glm::vec3& bar, glm::vec4& color) override
		{
			ZoneScoped;

			// Interpolate UV
			glm::vec2 uv = bar.x * varying_uv[0] +
				bar.y * varying_uv[1] +
				bar.z * varying_uv[2];

			// Normal map
			glm::vec4 normMapValue = material->bumpTexture->sample(uv);
			glm::vec3 normalMapSample = glm::vec3(normMapValue) * 2.0f - glm::vec3(1.0f);
			normalMapSample = glm::normalize(normalMapSample);

			// Interpolate TBN
			glm::mat3 interpolatedTBN =
				bar.x * varying_tbn0
				+ bar.y * varying_tbn1
				+ bar.z * varying_tbn2;

			// Extract T, B, N
			glm::vec3 interpolatedT = interpolatedTBN[0];
			glm::vec3 interpolatedB = interpolatedTBN[1];
			glm::vec3 interpolatedN = interpolatedTBN[2];

			// Re-orthogonalize
			interpolatedN = glm::normalize(interpolatedN);

			glm::vec3 interpolatedTNorm = glm::normalize(
				interpolatedT - interpolatedN * glm::dot(interpolatedN, interpolatedT)
			);
			glm::vec3 interpolatedBNorm = glm::normalize(
				interpolatedB
				- glm::dot(interpolatedB, interpolatedT) * interpolatedT
				- glm::dot(interpolatedB, interpolatedN) * interpolatedN
			);

			// Recompute final N from T x B if desired:
			// interpolatedN = glm::normalize(glm::cross(interpolatedTNorm, interpolatedB));

			// Reconstruct TBN
			glm::mat3 finalInterpolatedTBN;
			finalInterpolatedTBN[0] = interpolatedTNorm;
			finalInterpolatedTBN[1] = interpolatedBNorm;
			finalInterpolatedTBN[2] = interpolatedN;

			// Textures
			glm::vec3 albedo = glm::vec3(material->diffuseTexture->sample(uv));
			// Convert from sRGB to linear space
			albedo = pow(albedo, glm::vec3(2.2f));

			glm::vec4 metallicVec = material->metallicTexture->sample(uv);
			glm::vec4 roughnessVec = material->roughnessTexture ?
				material->roughnessTexture->sample(uv) :
				glm::vec4(1.0f);
			glm::vec4 aoVec = material->aoTexture ?
				material->aoTexture->sample(uv) :
				glm::vec4(1.0f);

			float metallic = metallicVec.x;
			float roughness = roughnessVec.x;
			float ao = aoVec.x;

			// Final normal in world space
			glm::vec3 normal = glm::normalize(finalInterpolatedTBN * normalMapSample);

			// Lighting setup
			glm::vec3 fragPos = bar.x * varying_fragWorldPos[0] +
				bar.y * varying_fragWorldPos[1] +
				bar.z * varying_fragWorldPos[2];

			glm::vec3 viewDir = glm::normalize(cameraPosition - fragPos);
			glm::vec3 lightDir = -lightDirection;
			glm::vec3 halfwayDir = glm::normalize(lightDir + viewDir);

			// Fresnel (Schlick)
			glm::vec3 F0 = lerp(glm::vec3(0.04f), albedo, metallic);
			glm::vec3 F = fresnelSchlick(std::max(glm::dot(halfwayDir, normal), 0.0f), F0);

			// NDF
			float NDF = distributionGGX(normal, halfwayDir, roughness);
			// Geometry
			float G = geometrySmith(normal, viewDir, lightDir, roughness);

			// Cook-Torrance
			glm::vec3 numerator = F * NDF * G;
			float denom = 4.0f * std::max(glm::dot(normal, viewDir), 0.0f)
				* std::max(glm::dot(normal, lightDir), 0.0f)
				+ 0.0001f;
			glm::vec3 specular = numerator / denom;

			// kS = F
			glm::vec3 kS = F;
			// kD = 1 - kS
			glm::vec3 kD = glm::vec3(1.0f) - kS;
			kD *= (1.0f - metallic);

			// Diffuse + specular
			float NdotL = std::max(glm::dot(normal, lightDir), 0.0f);
			glm::vec3 radiance = lightColor;
			glm::vec3 Lo = (kD * albedo / glm::pi<float>() + specular) * radiance * NdotL;

			// Ambient
			glm::vec3 ambient = glm::vec3(0.03f) * albedo * ao;

			// Combine
			glm::vec3 finalColor = ambient + Lo;

			// Tone-mapping (Reinhard)
			finalColor = finalColor / (finalColor + glm::vec3(1.0f));

			// Gamma correction
			finalColor = glm::pow(finalColor, glm::vec3(1.0f / 2.2f));
			color = glm::vec4(finalColor, 1.0f);

			return false;
		}

	public:
		glm::vec3 lightDirection;
		glm::vec3 lightColor;

	private:
		// --------------------------------------------------
		// Fresnel (Schlick)
		glm::vec3 fresnelSchlick(float cosTheta, const glm::vec3& F0)
		{
			return F0 + (glm::vec3(1.0f) - F0)
				* std::pow(std::clamp(1.0f - cosTheta, 0.0f, 1.0f), 5.0f);
		}

		// DistributionGGX (Trowbridge-Reitz)
		float distributionGGX(const glm::vec3& N, const glm::vec3& H, float roughness)
		{
			float a = roughness * roughness;
			float a2 = a * a;
			float NdotH = std::max(glm::dot(N, H), 0.0f);
			float NdotH2 = NdotH * NdotH;

			float num = a2;
			float denom = (NdotH2 * (a2 - 1.0f) + 1.0f);
			denom = glm::pi<float>() * denom * denom;
			return num / denom;
		}

		// Geometry Schlick-GGX
		float geometrySchlickGGX(float NdotV, float roughness)
		{
			float r = (roughness + 1.0f);
			float k = (r * r) / 8.0f;
			float num = NdotV;
			float denom = NdotV * (1.0f - k) + k;
			return num / denom;
		}

		// Geometry Smith
		float geometrySmith(const glm::vec3& N,
			const glm::vec3& V,
			const glm::vec3& L,
			float roughness)
		{
			float NdotV = std::max(glm::dot(N, V), 0.0f);
			float NdotL = std::max(glm::dot(N, L), 0.0f);
			float ggx1 = geometrySchlickGGX(NdotV, roughness);
			float ggx2 = geometrySchlickGGX(NdotL, roughness);
			return ggx1 * ggx2;
		}

		// Simple lerp
		glm::vec3 lerp(const glm::vec3& a, const glm::vec3& b, float t)
		{
			return a + t * (b - a);
		}

	private:
		glm::vec2 varying_uv[3];
		glm::vec3 varying_fragWorldPos[3];

		// TBN for each vertex of the triangle
		glm::mat3 varying_tbn0;
		glm::mat3 varying_tbn1;
		glm::mat3 varying_tbn2;
	};

} // namespace AR
