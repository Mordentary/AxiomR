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

			// Calculate MVP and Vertex Position
			glm::vec4 gl_Vertex = mvp * glm::vec4(vertex.position, 1.0f);

			// Store Varying UV and World Position - No changes needed
			varying_uv[indexInsideFace] = vertex.uv;
			varying_fragWorldPos[indexInsideFace] = glm::vec3(model * glm::vec4(vertex.position, 1.0f));

			// --- TBN Matrix Calculation ---

			// Optimize: Calculate TBN matrix elements directly in world space
			glm::vec3 T = glm::normalize(glm::vec3(model * glm::vec4(vertex.tangent, 0.0f)));
			glm::vec3 B = glm::normalize(glm::vec3(model * glm::vec4(vertex.bitangent, 0.0f)));
			glm::vec3 N = glm::normalize(glm::vec3(model * glm::vec4(vertex.normal, 0.0f)));

			glm::mat3 worldTBN;
			worldTBN[0] = T;
			worldTBN[1] = B;
			worldTBN[2] = N;

			if (indexInsideFace == 0) varying_tbn0 = worldTBN;
			else if (indexInsideFace == 1) varying_tbn1 = worldTBN;
			else if (indexInsideFace == 2) varying_tbn2 = worldTBN;

			return gl_Vertex;
		}

		virtual bool fragment(glm::vec3& bar, glm::vec4& color) override
		{
			ZoneScoped;

			glm::vec2 uv = bar.x * varying_uv[0] + bar.y * varying_uv[1] + bar.z * varying_uv[2];

			// --- Normal Mapping ---

			glm::vec4 normMapValue = material->bumpTexture->sample(uv);

			glm::vec3 normalMapSample = glm::normalize((glm::vec3(normMapValue) * 2.0f) - 1.0f);

			// --- TBN Interpolation ---

			glm::mat3 interpolatedTBN = bar.x * varying_tbn0 + bar.y * varying_tbn1 + bar.z * varying_tbn2;

			glm::vec3 interpolatedT = interpolatedTBN[0];
			glm::vec3 interpolatedB = interpolatedTBN[1];
			glm::vec3 interpolatedN = interpolatedTBN[2];

			// --- Re-orthogonalization ---

			interpolatedN = glm::normalize(interpolatedN); 
			glm::vec3 nt = interpolatedN * glm::dot(interpolatedN, interpolatedT);
			glm::vec3 interpolatedTNorm = glm::normalize(interpolatedT - nt);

			glm::vec3 interpolatedBNorm = glm::normalize(glm::cross(interpolatedN, interpolatedTNorm));

			interpolatedN = glm::normalize(glm::cross(interpolatedTNorm, interpolatedBNorm));

			// --- Texture Sampling and Pre-calculations ---

			glm::vec3 albedo = glm::vec3(material->diffuseTexture->sample(uv));
			albedo = pow(albedo, glm::vec3(2.2f)); // sRGB to linear

			// Optimize: Use conditional assignment (no branching) if textures might be null
			glm::vec3 metallicRoughnessAO =
				glm::vec3(material->metallicTexture->sample(uv).r, material->roughnessTexture->sample(uv).r, material->aoTexture->sample(uv).r);

			float metallic = metallicRoughnessAO.r;
			float roughness = metallicRoughnessAO.g;
			float ao = metallicRoughnessAO.b;

			float roughness2 = roughness * roughness;
			float roughness4 = roughness2 * roughness2;
			float oneMinusMetallic = 1.0f - metallic;

			// --- Lighting Calculations ---

			glm::vec3 normal = glm::normalize(interpolatedTBN * normalMapSample);

			glm::vec3 fragPos = bar.x * varying_fragWorldPos[0] + bar.y * varying_fragWorldPos[1] + bar.z * varying_fragWorldPos[2];
			glm::vec3 viewDir = glm::normalize(cameraPosition - fragPos);
			glm::vec3 lightDir = -lightDirection;
			glm::vec3 halfwayDir = glm::normalize(lightDir + viewDir);

			float NdotH = std::max(glm::dot(normal, halfwayDir), 0.0f);
			float NdotV = std::max(glm::dot(normal, viewDir), 0.0f);
			float NdotL = std::max(glm::dot(normal, lightDir), 0.0f);

			// --- Fresnel (Schlick) ---

			glm::vec3 F0 = lerp(glm::vec3(0.04f), albedo, metallic);

			// --- NDF (GGX) ---

			float NdotH2 = NdotH * NdotH;
			float denomPart = (NdotH2 * (roughness4 - 1.0f) + 1.0f);
			float NDF = roughness4 / (glm::pi<float>() * denomPart * denomPart);

			// --- Geometry (Smith) ---

			float r = roughness + 1.0f;
			float k = (r * r) / 8.0f;
			float NdotV_k = NdotV * (1.0f - k) + k;
			float NdotL_k = NdotL * (1.0f - k) + k;
			float G = (NdotV / NdotV_k) * (NdotL / NdotL_k);

			// --- Cook-Torrance BRDF ---

			glm::vec3 F = fresnelSchlick(std::max(glm::dot(halfwayDir, normal), 0.0f), F0);
			glm::vec3 numerator = F * NDF * G;
			float denom = 4.0f * NdotV * NdotL + 0.0001f;
			glm::vec3 specular = numerator / denom;

			// --- Diffuse and Ambient ---

			glm::vec3 kS = F;
			glm::vec3 kD = (glm::vec3(1.0f) - kS) * oneMinusMetallic;
			glm::vec3 diffuse = kD * albedo * (1.0f / glm::pi<float>());
			glm::vec3 ambient = glm::vec3(0.03f) * albedo * ao;

			// --- Combine and Post-processing ---

			glm::vec3 finalColor = ambient + (diffuse + specular) * lightColor * NdotL;

			// --- Tone Mapping and Gamma Correction ---
			finalColor = finalColor / (finalColor + glm::vec3(1.0f));
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
		inline glm::vec3 fresnelSchlick(float cosTheta, const glm::vec3& F0)
		{
			float oneMinusCosTheta = 1.0f - cosTheta;
			float term = oneMinusCosTheta * oneMinusCosTheta;
			term *= term;
			term *= oneMinusCosTheta; // Approximation of (1 - cosTheta)^5

			return F0 + (glm::vec3(1.0f) - F0) * term;
		}


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
}