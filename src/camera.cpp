#include "camera.hpp"
#include <algorithm>
#include <cmath>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/quaternion.hpp>

namespace AR {
	namespace {
		constexpr float PITCH_LIMIT = 89.0f;
		constexpr float MIN_DELTA = 0.0001f;
		constexpr float MAX_DELTA = 0.5f;

		float clampFov(float fov) {
			return std::clamp(fov, 1.0f, 179.0f);
		}

		float clampAspectRatio(float ratio) {
			return std::max(ratio, 0.1f);
		}

		void sanitizeRotationDelta(float& deltaYaw, float& deltaPitch) {
			deltaYaw = std::clamp(deltaYaw, -MAX_DELTA, MAX_DELTA);
			deltaPitch = std::clamp(deltaPitch, -MAX_DELTA, MAX_DELTA);
			if (std::abs(deltaYaw) < MIN_DELTA) deltaYaw = 0.0f;
			if (std::abs(deltaPitch) < MIN_DELTA) deltaPitch = 0.0f;
		}
	}

	Camera::Camera(const glm::vec3& position, const glm::vec3& target, float fov, float aspectRatio)
		: m_Position(position)
		, m_Fov(clampFov(fov))
		, m_AspectRatio(clampAspectRatio(aspectRatio))
	{
		glm::vec3 forward = glm::vec3(0.0f, 0.0f, -1.0f);
		glm::vec3 direction = glm::normalize(target - m_Position);

		float dot = glm::dot(forward, direction);

		if (glm::abs(dot + 1.0f) < 0.000001f) {
			m_Orientation = glm::angleAxis(glm::pi<float>(), glm::vec3(0.0f, 1.0f, 0.0f));
		}
		else {
			glm::vec3 rotationAxis = glm::cross(forward, direction);
			if (glm::length2(rotationAxis) < 0.000001f) {
				m_Orientation = glm::quat(1.0f, 0.0f, 0.0f, 0.0f);
			}
			else {
				float rotationAngle = std::acos(dot);
				m_Orientation = glm::angleAxis(rotationAngle, glm::normalize(rotationAxis));
			}
		}
		m_Orientation = glm::normalize(m_Orientation);
		updateVectors();
		updateProjectionMatrix();
	}

	void Camera::updateVectors() {
		m_Forward = glm::normalize(m_Orientation * glm::vec3(0.0f, 0.0f, -1.0f));
		m_Right = glm::normalize(m_Orientation * glm::vec3(1.0f, 0.0f, 0.0f));
		m_Up = glm::normalize(m_Orientation * glm::vec3(0.0f, 1.0f, 0.0f));
		m_ViewDirty = true;
	}

	void Camera::updateViewMatrix() {
		if (m_ViewDirty) {
			glm::mat4 rotationMatrix = glm::mat4_cast(glm::inverse(m_Orientation));
			glm::mat4 translationMatrix = glm::translate(glm::mat4(1.0f), -m_Position);
			m_View = rotationMatrix * translationMatrix;
			m_ViewProjection = m_Projection * m_View;
			m_ViewDirty = false;
		}
	}

	void Camera::updateProjectionMatrix() {
		if (m_ProjectionDirty) {
			m_Projection = glm::perspective(glm::radians(m_Fov), m_AspectRatio, m_NearPlane, m_FarPlane);
			m_ViewProjection = m_Projection * m_View;
			m_ProjectionDirty = false;
		}
	}

	void Camera::update(float deltaTime)
	{
		if (m_DeltaPitch != 0.f && m_DeltaYaw != 0.f)
		{
			float deltaYaw = glm::radians(m_DeltaYaw * deltaTime);
			float deltaPitch = glm::radians(m_DeltaPitch * deltaTime);
			sanitizeRotationDelta(deltaYaw, deltaPitch);
			// Apply rotation with smoothing
			applyRotation(deltaYaw, deltaPitch);
		}

		if (glm::length2(m_CurrentVelocity) > 0.0f) {
			m_Position += m_CurrentVelocity * deltaTime;
		}
		updateViewMatrix();
	}

	void Camera::handleInput(const Window& window) {
		if (window.isMouseButtonDown(1)) { // Right mouse button
			glm::vec2 mousePos = window.getMousePos();
			m_DeltaYaw = -(mousePos.x - m_LastMousePos.x) * m_MouseSensitivity;
			m_DeltaPitch = (mousePos.y - m_LastMousePos.y) * m_MouseSensitivity;
			m_LastMousePos = mousePos;
			m_ViewDirty = true;
		}
		else {
			m_DeltaPitch = 0.f;
			m_DeltaYaw = 0.f;
			m_LastMousePos = window.getMousePos(); // Update even when not rotating
		}

		// Keyboard movement
		glm::vec3 velocity(0.0f);
		if (window.isKeyDown('W'))
			velocity += m_Forward;
		if (window.isKeyDown('S'))
			velocity -= m_Forward;
		if (window.isKeyDown('A'))
			velocity -= m_Right;
		if (window.isKeyDown('D'))
			velocity += m_Right;

		m_CurrentVelocity = velocity;
	}

	void Camera::applyRotation(float deltaYaw, float deltaPitch) {
		constexpr float maxRotationSpeed = glm::radians(45.0f);
		deltaYaw = glm::clamp(deltaYaw, -maxRotationSpeed, maxRotationSpeed);
		deltaPitch = glm::clamp(deltaPitch, -maxRotationSpeed, maxRotationSpeed);

		glm::quat yawQuat = glm::angleAxis(deltaYaw, glm::vec3(0.0f, 1.0f, 0.0f));
		glm::quat pitchQuat = glm::angleAxis(deltaPitch, m_Right);

		glm::quat newOrientation = yawQuat * pitchQuat * m_Orientation;

		glm::vec3 newUp = glm::normalize(newOrientation * glm::vec3(0.0f, 1.0f, 0.0f));
		float upDot = glm::dot(newUp, glm::vec3(0.0f, 1.0f, 0.0f));

		if (upDot < 0.0f) {
			m_Orientation = glm::normalize(yawQuat * m_Orientation);
		}
		else {
			glm::vec3 newForward = glm::normalize(newOrientation * glm::vec3(0.0f, 0.0f, -1.0f));
			float pitch = glm::degrees(glm::asin(glm::clamp(newForward.y, -1.0f, 1.0f)));

			// Apply rotation based on pitch limits
			if (std::abs(pitch) <= PITCH_LIMIT) {
				m_Orientation = glm::normalize(newOrientation);
			}
			else {
				// Only apply yaw if we would exceed pitch limits
				m_Orientation = glm::normalize(yawQuat * m_Orientation);
			}
		}

		updateVectors();
	}

	void Camera::setAspectRatio(float aspectRatio) {
		m_AspectRatio = clampAspectRatio(aspectRatio);
		m_ProjectionDirty = true;
	}

	void Camera::setPosition(const glm::vec3& position) {
		m_Position = position;
		m_ViewDirty = true;
	}

	void Camera::setTarget(const glm::vec3& target) {
		glm::vec3 forward = glm::vec3(0.0f, 0.0f, -1.0f);
		glm::vec3 direction = glm::normalize(target - m_Position);
		m_Orientation = glm::rotation(forward, direction);
		updateVectors();
	}

	void Camera::setFov(float fov) {
		m_Fov = clampFov(fov);
		m_ProjectionDirty = true;
	}
	void Camera::setViewport(int x, int y, int width, int height) {
		m_ViewportX = x;
		m_ViewportY = y;
		m_ViewportWidth = width;
		m_ViewportHeight = height;
	}
	glm::mat4 Camera::getViewportMatrix() const {
		glm::mat4 viewport = glm::mat4(1.0f);

		float halfWidth = m_ViewportWidth * 0.5f;
		float halfHeight = m_ViewportHeight * 0.5f;

		viewport[0][0] = halfWidth;
		viewport[1][1] = -halfHeight; // Flip Y-axis to match screen coordinates
		viewport[2][2] = 0.5f;
		viewport[3][0] = m_ViewportX + halfWidth;
		viewport[3][1] = m_ViewportY + halfHeight;
		viewport[3][2] = 0.5f;

		return viewport;
	}
}