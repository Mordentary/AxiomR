#include "camera.hpp"

#include"window.hpp"
#include <cmath>

namespace AR {
	// Constructor
	Camera::Camera(
		Vec3f position,
		Vec3f target,
		Vec3f up,
		float fovRadians,
		float aspect,
		int viewportWidth,
		int viewportHeight,
		int viewportX,
		int viewportY,
		float nearZ,
		float farZ)
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
		m_ViewportHeight(viewportHeight),
		m_Yaw(-1.57f),  // Initialize yaw to -90 degrees (facing -Z direction)
		m_Pitch(0.0f)
	{
		updateVectors(); // Calculate initial front, right, and up vectors
	}

	// Setters
	void Camera::setPosition(const Vec3f& pos) {
		m_Position = pos;
		updateVectors();
	}

	void Camera::setTarget(const Vec3f& tgt) {
		m_Target = tgt;
		updateVectors();
	}

	void Camera::setUp(const Vec3f& u) {
		m_Up = u;
		updateVectors();
	}

	void Camera::setPerspective(float fov, float a, float n, float f) {
		m_FovRadians = fov;
		m_Aspect = a;
		m_NearZ = n;
		m_FarZ = f;
	}

	void Camera::setViewport(int x, int y, int width, int height) {
		m_ViewportX = x;
		m_ViewportY = y;
		m_ViewportWidth = width;
		m_ViewportHeight = height;
	}

	// Camera Control
	void Camera::rotate(float yaw, float pitch) {
		m_Yaw += yaw;
		m_Pitch += pitch;

		// Clamp pitch to avoid flipping the camera over
		if (m_Pitch > 1.5f) {
			m_Pitch = 1.5f;
		}
		if (m_Pitch < -1.5f) {
			m_Pitch = -1.5f;
		}

		updateVectors();
	}

	void Camera::handleInput(const Window& window, float deltaTime) {
		// Rotation
		if (window.isMouseButtonDown(1)) { // 0 = right mouse button
			Vec2i mousePos = window.getMousePos();
			float deltaX = (mousePos.x - m_LastMousePos.x) * m_MouseSensitivity * deltaTime;
			float deltaY = (mousePos.y - m_LastMousePos.y) * m_MouseSensitivity * deltaTime;

			rotate(deltaX, deltaY);

			m_LastMousePos = mousePos;
		}
		else {
			m_LastMousePos = window.getMousePos();
		}

		if (window.isKeyDown('W')) {
			move(m_Front, m_MovementSpeed * deltaTime);
		}
		if (window.isKeyDown('S')) {
			move(m_Front, -m_MovementSpeed * deltaTime);
		}
		if (window.isKeyDown('A')) {
			move(m_Right, -m_MovementSpeed * deltaTime);
		}
		if (window.isKeyDown('D')) {
			move(m_Right, m_MovementSpeed * deltaTime);
		}
	}

	void Camera::move(const Vec3f& direction, float amount) {
		m_Position += direction * amount;
		m_Target += direction * amount;

		updateVectors();
	}

	// Matrix Generation
	mat4f Camera::getViewMatrix() const {
		return mat4f::lookAt(m_Position, m_Target, m_Up);
	}

	mat4f Camera::getProjectionMatrix() const {
		return mat4f::perspective(m_FovRadians, m_Aspect, m_NearZ, m_FarZ);
	}

	mat4f Camera::getViewportMatrix() const {
		mat4f viewport = mat4f::identity();

		float halfWidth = m_ViewportWidth / 2.0f;
		float halfHeight = m_ViewportHeight / 2.0f;

		viewport[0][0] = halfWidth;
		viewport[1][1] = -halfHeight; // Flip Y-axis to match screen coordinates
		viewport[2][2] = 0.5f;
		viewport[3][0] = m_ViewportX + halfWidth;
		viewport[3][1] = m_ViewportY + halfHeight;
		viewport[3][2] = 0.5f;

		return viewport;
	}

	// Private Helper Function
	void Camera::updateVectors() {
		Vec3f front;
		front.x = cos(m_Yaw) * cos(m_Pitch);
		front.y = sin(m_Pitch);
		front.z = sin(m_Yaw) * cos(m_Pitch);
		m_Front = front.normalized();

		m_Right = m_Front.cross(Vec3f{ 0, 1, 0 }).normalized();
		m_Up = m_Right.cross(m_Front).normalized();

		m_Target = m_Position + m_Front;
	}
}