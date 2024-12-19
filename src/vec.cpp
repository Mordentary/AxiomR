#include "vec.hpp"
namespace AR
{
	Vec3f barycentric(const std::array<Vec2i, 3>& pts, const Vec2i& P) {
		Vec3f u;
		float det = (pts[1].y - pts[2].y) * (pts[0].x - pts[2].x) +
			(pts[2].x - pts[1].x) * (pts[0].y - pts[2].y);
		if (det == 0.0) return Vec3f(-1, 1, 1); // Degenerate triangle

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
}