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

	struct DefaultShader : public IShader {
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
			Vec3f lightDir = varying_lightDir * bar;
			float intensity = lightDir.dot(normal);
			intensity = std::clamp(intensity, 0.0f, 1.0f);
			color = material->diffuseTexture->sample(uv) * intensity;
			return false;
		}
	public:
		Vec3f lightDirection;
	private:
		Mat<2, 3, float> varying_uv;
		Mat<3, 3, float> varying_lightDir;
	};
}