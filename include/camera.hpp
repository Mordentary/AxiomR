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
			float aspect = 16 / 9,
			float nearZ = 0.1f,
			float farZ = 100.0f)
			: m_Position(position),
			m_Target(target),
			m_Up(up),
			m_FovRadians(fovRadians),
			m_Aspect(aspect),
			m_NearZ(nearZ),
			m_FarZ(farZ) {
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
		Mat4f getViewMatrix() const {
			return Mat4f::lookAt(m_Position, m_Target, m_Up);
		}

		Mat4f getProjectionMatrix() const {
			return Mat4f::perspective(m_FovRadians, m_Aspect, m_NearZ, m_FarZ);
		}
	private:
		Vec3f m_Position;
		Vec3f m_Target;
		Vec3f m_Up;
		float m_FovRadians;
		float m_Aspect;
		float m_NearZ;
		float m_FarZ;
	};
}