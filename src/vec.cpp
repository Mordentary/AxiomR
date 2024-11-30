#include "vec.hpp"
namespace AR
{
	Vec3f barycentric(Vec2i* pts, Vec2i P) {
		Vec3f u = Vec3f(pts[2].x - pts[0].x, pts[1].x - pts[0].x, pts[0].x - P.x)
			.cross(Vec3f(pts[2].y - pts[0].y, pts[1].y - pts[0].y, pts[0].y - P.y));

		if (std::abs(u.z) < 1) return Vec3f(-1, 1, 1);

		return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
	}
}