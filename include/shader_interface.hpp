#pragma once
#include"pipeline.hpp"
#include"vec.hpp"

namespace AR
{
	class Pipeline;
	struct Vertex;
	struct IShader {
		friend Pipeline;
	public:
		virtual ~IShader() = default;
	private:
		virtual Vec4f vertex(const Vertex& vertex, int indexInsideFace) = 0;
		virtual bool fragment(const Vec3f& bar, Vec4f& color) = 0;
	protected:
		Mat<4, 4, float> mvp;
		Mat<4, 4, float> viewportMat;
		Pipeline* m_PipelineState;
	};
}