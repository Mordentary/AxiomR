#pragma once
#include"math.hpp"
#include"mesh.hpp"

namespace AR
{
	class Pipeline;
	class TiledPipeline;
	struct Vertex;

	struct VertexOutput {
		glm::vec3 normal;
		float zNDC = 1.0f;
		glm::vec2 uv;
		glm::vec3 worldPos;
		glm::mat3 tbn;
	};

	struct VSTransformedTriangle
	{
		VertexOutput vertices[3];
		VertexOutput& operator [](int idx) {
			return vertices[idx];
		}
		VertexOutput operator [](int idx) const {
			return vertices[idx];
		}
	};

	struct IShader {
		friend Pipeline;
		friend TiledPipeline;
	public:
		virtual ~IShader() = default;
	private:
		virtual VertexOutput vertex(const Vertex& vertex, int indexInsideFace) = 0;
		virtual bool fragment(glm::vec3& bar, glm::vec4& color, const VSTransformedTriangle& tri) = 0;
	protected:
		glm::mat4 mvp;
		glm::mat4 model;
		glm::mat4 viewProj;
		glm::mat4 viewportMat;
		glm::vec3 cameraPosition;
		Pipeline* pipelineState;
		const Material* material;
	};
}