#pragma once
#include"pipeline.hpp"
#include"vec.hpp"

namespace AR
{
	class Pipeline;
	struct Vertex;
	struct IShader {
		virtual ~IShader() = default;
		virtual Vec4i vertex(const Vertex& vertex, int indexInsideFace) = 0;
		virtual bool fragment(const Vec3f& bar, Vec4f& color) = 0;
	protected:
		Pipeline* m_PipelineState;
	};
}