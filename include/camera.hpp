#pragma once
#include "vec.hpp"

namespace AR {
	class Window;
	class Camera {
	public:
		Camera(
			Vec3f position = { 0.0f, 0.0f, 3.0f },
			Vec3f target = { 0.0f, 0.0f, 0.0f },
			Vec3f up = { 0.0f, 1.0f, 0.0f },
			float fovRadians = 1.0472f,
			float aspect = 16.0f / 9.0f,
			int viewportWidth = 800,
			int viewportHeight = 600,
			int viewportX = 0,
			int viewportY = 0,
			float nearZ = 0.1f,
			float farZ = 100.0f
		);

		void setPosition(const Vec3f& pos);
		void setTarget(const Vec3f& tgt);
		void setUp(const Vec3f& u);
		void setPerspective(float fov, float a, float n, float f);
		void setViewport(int x, int y, int width, int height);

		Vec3f getPosition() const { return m_Position; }
		Vec3f getTarget() const { return m_Target; }
		Vec3f getUp() const { return m_Up; }
		Vec3f getFront() const { return m_Front; }
		Vec3f getRight() const { return m_Right; }
		float getYaw() const { return m_Yaw; }
		float getPitch() const { return m_Pitch; }

		void handleInput(const Window& window, float deltaTime);

		mat4f getViewMatrix() const;
		mat4f getProjectionMatrix() const;
		mat4f getViewportMatrix() const;

	private:
		void updateVectors();
		void rotate(float yaw, float pitch);
		void move(const Vec3f& direction, float amount);
	private:
		// Camera Attributes
		Vec3f m_Position;
		Vec3f m_Target;
		Vec3f m_Up;
		Vec3f m_Front;
		Vec3f m_Right;
		float m_Yaw;
		float m_Pitch;
		float m_MovementSpeed = 4.0;

		// Projection Parameters
		float m_FovRadians;
		float m_Aspect;
		float m_NearZ;
		float m_FarZ;

		// Viewport Parameters
		int m_ViewportX;
		int m_ViewportY;
		int m_ViewportWidth;
		int m_ViewportHeight;

		// Input Handling
		Vec2i m_LastMousePos;
	};
}