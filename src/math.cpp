#include "math.hpp"
namespace AR
{
	glm::vec3 barycentric(const std::array<glm::uvec2, 3>& pts, const glm::uvec2& P) {
		glm::vec3 u;
		float det = (pts[1].y - pts[2].y) * (pts[0].x - pts[2].x) +
			(pts[2].x - pts[1].x) * (pts[0].y - pts[2].y);
		if (det == 0.0) return glm::vec3(-1, 1, 1); // Degenerate triangle

		u.x = ((pts[1].y - pts[2].y) * (P.x - pts[2].x) +
			(pts[2].x - pts[1].x) * (P.y - pts[2].y)) / det;
		u.y = ((pts[2].y - pts[0].y) * (P.x - pts[2].x) +
			(pts[0].x - pts[2].x) * (P.y - pts[2].y)) / det;
		u.z = 1.0f - u.x - u.y;
		return u;
	}

	double degreeToRad(double degrees)
	{
		return degrees * std::numbers::pi / 180.0;
	}

	glm::vec3 reflect(const glm::vec3& incident, const glm::vec3& normal) {
		return incident - 2.0f * normal * glm::dot(normal,incident);
	}
}