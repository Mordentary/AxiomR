#pragma once
#include <cmath>
#include <iostream>
#include <type_traits>

namespace AR {
	//============================================================
	// Vec2: 2D Vector
	//============================================================
	template<typename T>
	class Vec2 {
		static_assert(std::is_arithmetic<T>::value, "Vec2 requires an arithmetic type.");
	public:
		T x, y;

		// Constructors
		constexpr Vec2() : x(T(0)), y(T(0)) {}
		constexpr Vec2(T xVal, T yVal) : x(xVal), y(yVal) {}
		constexpr Vec2(const Vec2& v) = default;
		Vec2& operator=(const Vec2& v) = default;

		// Dot product
		constexpr T dot(const Vec2& v) const {
			return x * v.x + y * v.y;
		}

		// Squared length
		constexpr T lengthSquared() const {
			return x * x + y * y;
		}

		// Length (magnitude)
		T length() const {
			return std::sqrt(lengthSquared());
		}

		// Normalized vector (returns a new vector)
		Vec2 normalized() const {
			T len = length();
			if (len > T(0)) {
				return Vec2(x / len, y / len);
			}
			// If length is zero, return the same vector (or consider returning a default unit vector)
			return *this;
		}

		// Normalize this vector in-place
		void normalize() {
			T len = length();
			if (len > T(0)) {
				x /= len;
				y /= len;
			}
		}

		// Arithmetic operators
		constexpr Vec2 operator+(const Vec2& v) const {
			return Vec2(x + v.x, y + v.y);
		}

		constexpr Vec2 operator-(const Vec2& v) const {
			return Vec2(x - v.x, y - v.y);
		}

		constexpr Vec2 operator*(T scalar) const {
			return Vec2(x * scalar, y * scalar);
		}

		constexpr Vec2 operator/(T scalar) const {
			return Vec2(x / scalar, y / scalar);
		}

		constexpr Vec2& operator+=(const Vec2& v) {
			x += v.x;
			y += v.y;
			return *this;
		}

		constexpr Vec2& operator-=(const Vec2& v) {
			x -= v.x;
			y -= v.y;
			return *this;
		}

		constexpr Vec2& operator*=(T scalar) {
			x *= scalar;
			y *= scalar;
			return *this;
		}

		constexpr Vec2& operator/=(T scalar) {
			x /= scalar;
			y /= scalar;
			return *this;
		}

		constexpr bool operator==(const Vec2& v) const {
			return (x == v.x) && (y == v.y);
		}

		constexpr bool operator!=(const Vec2& v) const {
			return !(*this == v);
		}
	};

	// Scalar multiplication (scalar * vector)
	template<typename T>
	constexpr Vec2<T> operator*(T scalar, const Vec2<T>& v) {
		return v * scalar;
	}

	//============================================================
	// Vec3: 3D Vector
	//============================================================
	template<typename T>
	class Vec3 {
		static_assert(std::is_arithmetic<T>::value, "Vec3 requires an arithmetic type.");
	public:
		T x, y, z;

		constexpr Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
		constexpr Vec3(T xVal, T yVal, T zVal) : x(xVal), y(yVal), z(zVal) {}
		constexpr Vec3(const Vec3& v) = default;
		Vec3& operator=(const Vec3& v) = default;

		// Dot product
		constexpr T dot(const Vec3& v) const {
			return x * v.x + y * v.y + z * v.z;
		}

		// Cross product
		constexpr Vec3 cross(const Vec3& v) const {
			return Vec3(
				(y * v.z - z * v.y),
				(z * v.x - x * v.z),
				(x * v.y - y * v.x)
			);
		}

		// Squared length
		constexpr T lengthSquared() const {
			return x * x + y * y + z * z;
		}

		// Length (magnitude)
		T length() const {
			return std::sqrt(lengthSquared());
		}

		// Normalized vector (returns a new vector)
		Vec3 normalized() const {
			T len = length();
			if (len > T(0)) {
				return Vec3(x / len, y / len, z / len);
			}
			// If length is zero, return the same vector (or consider returning a default unit vector)
			return *this;
		}

		// Normalize this vector in-place
		void normalize() {
			T len = length();
			if (len > T(0)) {
				x /= len;
				y /= len;
				z /= len;
			}
		}

		// Arithmetic operators
		constexpr Vec3 operator+(const Vec3& v) const {
			return Vec3(x + v.x, y + v.y, z + v.z);
		}

		constexpr Vec3 operator-(const Vec3& v) const {
			return Vec3(x - v.x, y - v.y, z - v.z);
		}

		constexpr Vec3 operator*(T scalar) const {
			return Vec3(x * scalar, y * scalar, z * scalar);
		}

		constexpr Vec3 operator/(T scalar) const {
			return Vec3(x / scalar, y / scalar, z / scalar);
		}

		constexpr Vec3 operator-() const {
			return Vec3(-x, -y, -z);
		}

		constexpr Vec3& operator+=(const Vec3& v) {
			x += v.x;
			y += v.y;
			z += v.z;
			return *this;
		}

		constexpr Vec3& operator-=(const Vec3& v) {
			x -= v.x;
			y -= v.y;
			z -= v.z;
			return *this;
		}

		constexpr Vec3& operator*=(T scalar) {
			x *= scalar;
			y *= scalar;
			z *= scalar;
			return *this;
		}

		constexpr Vec3& operator/=(T scalar) {
			x /= scalar;
			y /= scalar;
			z /= scalar;
			return *this;
		}

		constexpr bool operator==(const Vec3& v) const {
			return (x == v.x) && (y == v.y) && (z == v.z);
		}

		constexpr bool operator!=(const Vec3& v) const {
			return !(*this == v);
		}
	};

	// Scalar multiplication (scalar * vector)
	template<typename T>
	constexpr Vec3<T> operator*(T scalar, const Vec3<T>& v) {
		return v * scalar;
	}

	//============================================================
	// Type Aliases
	//============================================================
	using Vec2f = Vec2<float>;
	using Vec2d = Vec2<double>;
	using Vec2i = Vec2<int>;

	using Vec3f = Vec3<float>;
	using Vec3d = Vec3<double>;
	using Vec3i = Vec3<int>;

	//============================================================
	// I/O Stream Operators
	//============================================================
	template<typename T>
	std::ostream& operator<<(std::ostream& os, const Vec2<T>& v) {
		return os << "Vec2(" << v.x << ", " << v.y << ")";
	}

	template<typename T>
	std::ostream& operator<<(std::ostream& os, const Vec3<T>& v) {
		return os << "Vec3(" << v.x << ", " << v.y << ", " << v.z << ")";
	}

	//============================================================
	// Forward Declaration of barycentric
	//============================================================
	struct Vec4f {
		float x, y, z, w;
	};

	Vec3f barycentric(Vec2i* pts, Vec2i P);

	//============================================================
	// Mat4f: 4x4 Matrix
	//============================================================
	struct Mat4f {
		float m[4][4];

		static Mat4f identity() {
			Mat4f I;
			for (int r = 0; r < 4; r++)
				for (int c = 0; c < 4; c++)
					I.m[r][c] = (r == c) ? 1.0f : 0.0f;
			return I;
		}

		static Mat4f translate(const Vec3f& t) {
			Mat4f M = identity();
			M.m[0][3] = t.x;
			M.m[1][3] = t.y;
			M.m[2][3] = t.z;
			return M;
		}

		static Mat4f rotateX(float angleRadians) {
			Mat4f R = identity();
			float c = std::cos(angleRadians);
			float s = std::sin(angleRadians);
			R.m[1][1] = c;
			R.m[1][2] = -s;
			R.m[2][1] = s;
			R.m[2][2] = c;
			return R;
		}

		static Mat4f rotateY(float angleRadians) {
			Mat4f R = identity();
			float c = std::cos(angleRadians);
			float s = std::sin(angleRadians);
			R.m[0][0] = c;
			R.m[0][2] = s;
			R.m[2][0] = -s;
			R.m[2][2] = c;
			return R;
		}

		static Mat4f rotateZ(float angleRadians) {
			Mat4f R = identity();
			float c = std::cos(angleRadians);
			float s = std::sin(angleRadians);
			R.m[0][0] = c;
			R.m[0][1] = -s;
			R.m[1][0] = s;
			R.m[1][1] = c;
			return R;
		}

		static Mat4f rotateXYZ(float angleX, float angleY, float angleZ) {
			Mat4f rx = rotateX(angleX);
			Mat4f ry = rotateY(angleY);
			Mat4f rz = rotateZ(angleZ);

			return rz * ry * rx;
		}

		static Mat4f lookAt(const Vec3f& eye, const Vec3f& target, const Vec3f& up) {
			Vec3f f = target - eye;
			f.normalize();
			Vec3f r = f.cross(up);
			r.normalize();
			Vec3f u = -r.cross(f);

			Mat4f view = identity();
			view.m[0][0] = r.x; view.m[0][1] = r.y; view.m[0][2] = r.z; view.m[0][3] = -r.dot(eye);
			view.m[1][0] = u.x; view.m[1][1] = u.y; view.m[1][2] = u.z; view.m[1][3] = -u.dot(eye);
			view.m[2][0] = f.x; view.m[2][1] = f.y; view.m[2][2] = f.z; view.m[2][3] = -f.dot(eye);
			view.m[3][0] = 0.0f; view.m[3][1] = 0.0f; view.m[3][2] = 0.0f; view.m[3][3] = 1.0f;
			return view;
		}

		static Mat4f perspective(float fovRadians, float aspect, float nearZ, float farZ) {
			float f = 1.0f / std::tan(fovRadians / 2.0f);
			Mat4f P = {};
			P.m[0][0] = f / aspect;
			P.m[1][1] = f;
			P.m[2][2] = (farZ + nearZ) / (nearZ - farZ);
			P.m[2][3] = (2 * farZ * nearZ) / (nearZ - farZ);
			P.m[3][2] = -1.0f;
			return P;
		}

		Vec4f operator*(const Vec4f& v) const {
			Vec4f res;
			res.x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3] * v.w;
			res.y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3] * v.w;
			res.z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3] * v.w;
			res.w = m[3][0] * v.x + m[3][1] * v.y + m[3][2] * v.z + m[3][3] * v.w;
			return res;
		}

		Mat4f operator*(const Mat4f& other) const {
			Mat4f result = {};
			for (int r = 0; r < 4; r++) {
				for (int c = 0; c < 4; c++) {
					float sum = 0.0f;
					for (int k = 0; k < 4; k++) {
						sum += m[r][k] * other.m[k][c];
					}
					result.m[r][c] = sum;
				}
			}
			return result;
		}
	};
}