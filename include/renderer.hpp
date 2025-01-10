#pragma once
#include <vector>
#include "camera.hpp"
#include"mesh.hpp"
#include"window.hpp"
#include"IShader.hpp"
#include"pipeline.hpp"
#include"timer.hpp"
#include"glm\glm.hpp"
namespace AR
{
	class Vertex;
	class Pipeline;

	class Renderer
	{
	public:
		explicit Renderer() = default;
		~Renderer();
		void init(uint32_t screenWidth, uint32_t screenHeight);
		void run();

		void drawLine(glm::vec2 p1, glm::vec2 p2, Color color);
		void drawTriangle(const Vertex& p0, const  Vertex& p1, const  Vertex& p2, glm::vec3 normal, glm::vec3 ligthDir);
		void drawMesh(const glm::mat4& trnsfrm, const Mesh& mesh);
	private:
		void render();
	private:
		std::unique_ptr<Framebuffer> m_Framebuffer;
		std::unique_ptr<Window> m_Window;
		std::unique_ptr<Camera> m_Camera;
		std::unique_ptr<Pipeline> m_DefaultPipeline;
		std::vector<std::unique_ptr<Mesh>> m_Meshes;
		std::vector<std::unique_ptr<IShader>> m_Shaders;
		IShader* m_CurrentShader;
		Timer m_FrameTime;
	};
}