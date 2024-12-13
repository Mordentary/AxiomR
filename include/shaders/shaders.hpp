#include"../pipeline.hpp"

namespace AR
{
	struct FlatShader : public IShader {
		friend Pipeline;
	public:
		virtual ~FlatShader() {}

		virtual Vec4f vertex(const Vertex& vertex, int indexInsideFace) {
			Vec4f gl_Vertex = toVec4f(vertex.position, 1.0f);
			gl_Vertex = mvp * gl_Vertex;
			interpolated_triangle.set_col(indexInsideFace, perspectiveDivision(gl_Vertex));
			return (gl_Vertex);
		}

		virtual bool fragment(const Vec3f& bar, Vec4f& color)
		{
			Vec3f n = (interpolated_triangle.get_col(1) - interpolated_triangle.get_col(0)).cross(interpolated_triangle.get_col(2) - interpolated_triangle.get_col(0));
			n.normalize();
			float intensity = lightDirection.dot(n);
			intensity = std::clamp(intensity, 0.0f, 1.0f);
			color = Vec4f{ 1.0f, 1.0f, 1.0f, 1.0f } *intensity;;
			return false;
		}

	public:
		Vec3f lightDirection;
	private:
		Mat<3, 3, float> interpolated_triangle;
	};
}
