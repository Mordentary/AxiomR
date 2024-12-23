#include "pipeline.hpp"
#include "IShader.hpp"
#include "camera.hpp"
#include "framebuffer.hpp"
namespace AR
{
	Pipeline::Pipeline()
		: m_Shader(nullptr), m_Camera(nullptr), m_Framebuffer(nullptr)
	{
	}

	void Pipeline::setShader(IShader* shader) {
		m_Shader = shader;
	}

	void Pipeline::setCamera(const Camera* cam) {
		m_Camera = cam;
	}

	void Pipeline::setFramebuffer(Framebuffer* fb) {
		m_Framebuffer = fb;
	}
	mat4f Pipeline::getViewportMat()
	{
		return m_Camera->getViewportMatrix();
	}

	void Pipeline::drawMesh(const mat4f& modelMatrix, const Mesh& mesh) {
		if (!m_Shader || !m_Camera || !m_Framebuffer) return;

		// Compute final transformation matrices once per draw
		mat4f view = m_Camera->getViewMatrix();
		mat4f proj = m_Camera->getProjectionMatrix();
		mat4f mvp = proj * view * modelMatrix;
		m_Shader->model = modelMatrix;
		m_Shader->viewProj = proj * view;
		m_Shader->mvp = m_Shader->viewProj * modelMatrix;
		m_Shader->viewportMat = m_Camera->getViewportMatrix();
		m_Shader->cameraPosition = m_Camera->getPosition();
		// Draw all faces
		const Vertex* vertices = mesh.getVertices().data();
		const std::vector<Face>& faces = mesh.getFaces();
		const auto& groups = mesh.getMaterialGroups();
		for (const auto& group : groups) {
			const Material* material = mesh.getMaterial(group.materialName);
			//renderer.setMaterial(material); // Set the material properties
			m_Shader->material = material;
			for (size_t i = group.startIndex; i < group.startIndex + group.faceCount; ++i) {
				const auto& face = faces[i];
				assert(face.vertexIndices.size() == 3);

				size_t v0 = face.vertexIndices[0];
				size_t v1 = face.vertexIndices[1];
				size_t v2 = face.vertexIndices[2];
				const Vertex& v0Data = vertices[v0];
				const Vertex& v1Data = vertices[v1];
				const Vertex& v2Data = vertices[v2];

				Vec4f clipCoords[3];
				clipCoords[0] = m_Shader->vertex(v0Data, 0);
				clipCoords[1] = m_Shader->vertex(v1Data, 1);
				clipCoords[2] = m_Shader->vertex(v2Data, 2);

				rasterizeTriangle(clipCoords);
			}
		}
	}

	Vec3f Pipeline::barycentric(const Vec2f& A, const Vec2f& B, const Vec2f& C, const Vec2f& P) const {
		Vec3f s0(C.x - A.x, B.x - A.x, A.x - P.x);
		Vec3f s1(C.y - A.y, B.y - A.y, A.y - P.y);
		Vec3f u = s0.cross(s1);
		if (std::abs(u.z) > 1e-2)
			return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
		return Vec3f(-1, 1, 1);
	}

	void Pipeline::rasterizeTriangle(Vec4f clip[3]) {
		int width = m_Framebuffer->getWidth();
		int height = m_Framebuffer->getHeight();

		// Perform perspective divide to get NDC coordinates
		Vec3f ndc[3];
		for (int i = 0; i < 3; i++) {
			float w = clip[i].w != 0.0f ? clip[i].w : 1.0f;
			ndc[i] = Vec3f(clip[i].x / w, clip[i].y / w, clip[i].z / w);
		}

		// Convert from NDC [-1,1] to screen coordinates [0,width], [0,height]
		Vec2i screen[3];
		for (int i = 0; i < 3; i++) {
			screen[i].x = static_cast<int>((ndc[i].x + 1.0f) * 0.5f * width);
			screen[i].y = static_cast<int>((ndc[i].y + 1.0f) * 0.5f * height);
		}

		// Compute bounding box
		Vec2i bboxmin(width - 1, height - 1);
		Vec2i bboxmax(0, 0);
		for (int i = 0; i < 3; i++) {
			bboxmin.x = std::max(0, std::min(bboxmin.x, screen[i].x));
			bboxmin.y = std::max(0, std::min(bboxmin.y, screen[i].y));
			bboxmax.x = std::min(width - 1, std::max(bboxmax.x, screen[i].x));
			bboxmax.y = std::min(height - 1, std::max(bboxmax.y, screen[i].y));
		}

		// Raster loop
		Vec2i P;
		for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
			for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
				Vec3f bc_screen = barycentric(
					Vec2f(static_cast<float>(screen[0].x), static_cast<float>(screen[0].y)),
					Vec2f(static_cast<float>(screen[1].x), static_cast<float>(screen[1].y)),
					Vec2f(static_cast<float>(screen[2].x), static_cast<float>(screen[2].y)),
					Vec2f(static_cast<float>(P.x), static_cast<float>(P.y))
				);
				if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;

				// Interpolate depth
				float z = ndc[0].z * bc_screen.x + ndc[1].z * bc_screen.y + ndc[2].z * bc_screen.z;

				// Depth test
				if (z < m_Framebuffer->getDepth(P.x, P.y)) {
					// Fragment shader
					Vec4f colorOut;
					bool discard = m_Shader->fragment(bc_screen, colorOut);
					if (!discard) {
						Color finalColor = {
							static_cast<uint8_t>(std::min(std::max(colorOut.x, 0.0f), 1.0f) * 255),
							static_cast<uint8_t>(std::min(std::max(colorOut.y, 0.0f), 1.0f) * 255),
							static_cast<uint8_t>(std::min(std::max(colorOut.z, 0.0f), 1.0f) * 255),
							static_cast<uint8_t>(std::min(std::max(colorOut.w, 0.0f), 1.0f) * 255)
						};
						m_Framebuffer->setDepth(P.x, P.y, z);
						m_Framebuffer->setPixel(P.x, P.y, finalColor);
					}
				}
			}
		}
	}
}