#pragma once
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include "window.hpp"
namespace AR {
	class Camera {
	public:
		Camera(const glm::vec3& position = { 0.0f, 0.0f, 0.0f },
			const glm::vec3& target = { 0.0f, 0.0f, -1.0f },
			float fov = 60.0f,
			float aspectRatio = 16.0f / 9.0f);

		void setPosition(const glm::vec3& position);
		void setTarget(const glm::vec3& target);
		void setFov(float fov);
		void setAspectRatio(float aspectRatio);
		void setViewport(int x, int y, int width, int height);
		// Accessors
		const glm::mat4& getViewMatrix() const { return m_View; }
		const glm::mat4& getProjectionMatrix() const { return m_Projection; }
		const glm::mat4& getViewProjectionMatrix() const { return m_ViewProjection; }
		glm::mat4 getViewportMatrix() const;
		glm::vec3 getPosition() const { return m_Position; };

		void update(float deltaTime);
		void handleInput(const Window& window);

	private:
		void updateVectors();
		void updateViewMatrix();
		void updateProjectionMatrix();
		void applyRotation(float deltaYaw, float deltaPitch);

	private:
		glm::vec3 m_Position;
		glm::quat m_Orientation;
		glm::vec3 m_Forward;
		glm::vec3 m_Right;
		glm::vec3 m_Up;

		float m_DeltaYaw, m_DeltaPitch;
		glm::vec3 m_CurrentVelocity;

		float m_Fov;
		float m_AspectRatio;
		float m_NearPlane = 0.1f;
		float m_FarPlane = 100.0f;

		float m_MovementSpeed = 5.0f;
		float m_MouseSensitivity = 1.0f;

		// For mouse input
		glm::vec2 m_LastMousePos = glm::vec2(0.0f);

		// Viewport parameters
		int m_ViewportX = 0;
		int m_ViewportY = 0;
		int m_ViewportWidth = 800;
		int m_ViewportHeight = 600;

		glm::mat4 m_View = glm::mat4(1.0f);
		glm::mat4 m_Projection = glm::mat4(1.0f);
		glm::mat4 m_ViewProjection = glm::mat4(1.0f);
		bool m_ViewDirty = true;
		bool m_ProjectionDirty = true;
	};
}