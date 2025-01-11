#include <camera.hpp>
#include <cmath>
#include <framebuffer.hpp>
#include <memory>
#include <mesh.hpp>
#include <pipeline.hpp>
#include <tiled_pipeline.hpp>
#include <shaders/shaders.hpp>
#include <utility>
#include <math.hpp>
#include <vector>
#include <window.hpp>
#include "renderer.hpp"
#include"glm\gtc\matrix_transform.hpp"
namespace AR {
	void Renderer::drawMesh(const glm::mat4& transMat, const Mesh& mesh)
	{
		m_DefaultPipeline->drawMesh(transMat, mesh);
	}

	void Renderer::drawLine(glm::vec2 p0, glm::vec2 p1, Color color)
	{
		int x0 = p0.x, y0 = p0.y;
		int x1 = p1.x, y1 = p1.y;

		bool steep = false;
		if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
			std::swap(x0, y0);
			std::swap(x1, y1);
			steep = true;
		}
		if (x0 > x1) {
			std::swap(x0, x1);
			std::swap(y0, y1);
		}

		int dx = x1 - x0;
		int dy = y1 - y0;
		int derror2 = std::abs(dy) * 2;
		int error2 = 0;
		int y = y0;

		for (int x = x0; x <= x1; x++) {
			if (steep) {
				m_Framebuffer.get()->setPixel(y, x, color);
			}
			else {
				m_Framebuffer.get()->setPixel(x, y, color);
			}
			error2 += derror2;
			if (error2 > dx) {
				y += (y1 > y0 ? 1 : -1);
				error2 -= dx * 2;
			}
		}
	}

	void Renderer::init(uint32_t screenWidth, uint32_t screenHeight)
	{
		//ZoneScoped;

		m_Window = std::make_unique<Window>(screenWidth, screenHeight, "AxiomR");
		m_Window->show();
		uint32_t renderAreaWidth = m_Window->m_Bitmap->getWidth();
		uint32_t renderAreaHeight = m_Window->m_Bitmap->getHeight();
		m_Framebuffer = std::make_unique<Framebuffer>(renderAreaWidth, renderAreaHeight, true);

		m_Camera = std::make_unique<Camera>(glm::vec3{ 0.0,0.0,5.0f }, glm::vec3{ 0.0,0.0,0.0f }, (60.f), (float)renderAreaWidth / renderAreaHeight);
		m_Camera->setViewport(0, 0, renderAreaWidth, renderAreaHeight);

		m_DefaultPipeline.reset(new TiledPipeline(std::thread::hardware_concurrency(), m_Camera.get(), m_Framebuffer.get()));

		FlatShader* flatShader = new FlatShader();
		flatShader->lightDirection = glm::vec3(0.0f, -1.0f, 0.0f);
		m_Shaders.emplace_back(flatShader);

		//DefaultShader* defaultShader = new DefaultShader();
		//defaultShader->lightDirection = glm::vec3(0.0f, -1.0f, 0.0f);
		//m_Shaders.emplace_back(defaultShader);

		PhongShader* phongShader = new PhongShader();
		phongShader->lightDirection = glm::vec3(0.0f, -1.0f, 0.0f);
		phongShader->lightColor = glm::vec3(0.6f, 0.6f, 0.6f);
		m_Shaders.emplace_back(phongShader);

		PBRShader* pbrShader = new PBRShader();
		pbrShader->lightDirection = glm::vec3(0.0f, -1.0f, 0.0f);
		pbrShader->lightColor = glm::vec3(5.6f, 5.6f, 5.6f);
		m_Shaders.emplace_back(pbrShader);

		m_CurrentShader = m_Shaders[2].get();
		m_DefaultPipeline->setShader(m_CurrentShader);
		m_Meshes.emplace_back(std::make_unique<Mesh>("assets/Gas Tank/Gas Tank.obj"));
	}
	Renderer::~Renderer()
	{
		m_Framebuffer.reset();
		m_Window.reset();
	}
	void Renderer::render()
	{
		glm::mat4 translation = glm::translate(glm::mat4(1.0f), glm::vec3(0, -3, 3));
		for (auto& mesh : m_Meshes)
		{
			float time = GetTickCount64() / 1000.0f;
			float angle = 0.4f * time;
			glm::mat4 rotation = glm::rotate(glm::mat4(1.0f), (angle), glm::vec3(0, 1, 0));
			drawMesh(translation * rotation, *mesh);
		}
	}


	void Renderer::run() {
		Color BLACK{ 0,0,0,255 };
		AR::DeltaTimeAverager deltaAverager(1.0f, "Averaged Delta Time");

		while (m_Window->isRunning()) {
			m_FrameTime.start();
			float deltaTimeS= m_FrameTime.getTimeSeconds();

			// Process input, update, render, etc.
			m_Window->processMessages();
			m_Camera->handleInput(*m_Window);
			m_Camera->update(deltaTimeS);
			m_Framebuffer->clearColor(BLACK);
			m_Framebuffer->clearDepth();
			render();
			m_Window->present(*m_Framebuffer);
			m_FrameTime.stop();
			deltaAverager.addDelta(deltaTimeS);

			FrameMark;
		}
	}
}