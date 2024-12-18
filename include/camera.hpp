#pragma once
#include"vec.hpp"

namespace AR
{
	class Camera
	{
	public:
		Camera(
			Vec3f position = { 0.0f, 0.0f, 3.0f },
			Vec3f target = { 0.0f, 0.0f, 0.0f },
			Vec3f up = { 0.0f, 1.0f, 0.0f },
			float fovRadians = 3.14159f / 4.0f,
			float aspect = 16.0f / 9.0f,  // Ensure floating point division
			float nearZ = 0.1f,
			float farZ = 100.0f,
			int viewportX = 0,
			int viewportY = 0,
			int viewportWidth = 800,
			int viewportHeight = 600)
			: m_Position(position),
			m_Target(target),
			m_Up(up),
			m_FovRadians(fovRadians),
			m_Aspect(aspect),
			m_NearZ(nearZ),
			m_FarZ(farZ),
			m_ViewportX(viewportX),
			m_ViewportY(viewportY),
			m_ViewportWidth(viewportWidth),
			m_ViewportHeight(viewportHeight)
		{
		}

		void setPosition(const Vec3f& pos) { m_Position = pos; }
		void setTarget(const Vec3f& tgt) { m_Target = tgt; }
		void setUp(const Vec3f& u) { m_Up = u; }
		void setPerspective(float fov, float a, float n, float f) {
			m_FovRadians = fov;
			m_Aspect = a;
			m_NearZ = n;
			m_FarZ = f;
		}
		void setViewport(int x, int y, int width, int height) {
			m_ViewportX = x;
			m_ViewportY = y;
			m_ViewportWidth = width;
			m_ViewportHeight = height;
		}
		mat4f getViewMatrix() const {
			return mat4f::lookAt(m_Position, m_Target, m_Up);
		}

		mat4f getProjectionMatrix() const {
			return mat4f::perspective(m_FovRadians, m_Aspect, m_NearZ, m_FarZ);
		}
		mat4f getViewportMatrix() const {
			mat4f viewport = mat4f::identity();

			float halfWidth = m_ViewportWidth / 2.0f;
			float halfHeight = m_ViewportHeight / 2.0f;

			viewport[0][0] = halfWidth;
			viewport[1][1] = halfHeight;
			viewport[2][2] = 0.5f;
			viewport[3][0] = m_ViewportX + halfWidth;
			viewport[3][1] = m_ViewportY + halfHeight;
			viewport[3][2] = 0.5f;

			return viewport;
		}
	private:
		Vec3f m_Position;
		Vec3f m_Target;
		Vec3f m_Up;
		float m_FovRadians;
		float m_Aspect;
		float m_NearZ;
		float m_FarZ;

		int m_ViewportX;
		int m_ViewportY;
		int m_ViewportWidth;
		int m_ViewportHeight;
	};
}