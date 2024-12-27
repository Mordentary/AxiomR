#pragma once
#include"math.hpp"
#include"mesh.hpp"

namespace AR
{
	class Pipeline;
	struct Vertex;

	struct IShader {
		friend Pipeline;
	public:
		virtual ~IShader() = default;
	private:
		virtual glm::vec4 vertex(const Vertex& vertex, int indexInsideFace) = 0;
		virtual bool fragment(glm::vec3& bar, glm::vec4& color) = 0;
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