#pragma once
#include <cmath>
#include <iostream>
#include <type_traits>

namespace AR {
	template<typename T>
	class Vec2 {
	public:
		T x, y;

		// Constructors
		constexpr Vec2() : x(0), y(0) {}
		constexpr Vec2(T x, T y) : x(x), y(y) {}
		constexpr Vec2(const Vec2& v) = default;

		// Assignment
		Vec2& operator=(const Vec2& v) = default;

		// Vector operations
		constexpr T dot(const Vec2& v) const {
			return x * v.x + y * v.y;
		}

		constexpr T lengthSquared() const {
			return x * x + y * y;
		}

		T length() const {
			return std::sqrt(lengthSquared());
		}

		Vec2 normalized() const {
			T len = length();
			if (len > 0) {
				return Vec2(x / len, y / len);
			}
			return *this;
		}

		void normalize() {
			T len = length();
			if (len > 0) {
				x /= len;
				y /= len;
			}
		}

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
			return x == v.x && y == v.y;
		}

		constexpr bool operator!=(const Vec2& v) const {
			return !(*this == v);
		}
	};

	template<typename T>
	class Vec3 {
	public:
		T x, y, z;

		constexpr Vec3() : x(0), y(0), z(0) {}
		constexpr Vec3(T x, T y, T z) : x(x), y(y), z(z) {}
		constexpr Vec3(const Vec3& v) = default;

		Vec3& operator=(const Vec3& v) = default;

		constexpr T dot(const Vec3& v) const {
			return x * v.x + y * v.y + z * v.z;
		}

		constexpr Vec3 cross(const Vec3& v) const {
			return Vec3(
				y * v.z - z * v.y,
				z * v.x - x * v.z,
				x * v.y - y * v.x
			);
		}

		constexpr T lengthSquared() const {
			return x * x + y * y + z * z;
		}

		T length() const {
			return std::sqrt(lengthSquared());
		}

		Vec3 normalized() const {
			T len = length();
			if (len > 0) {
				return Vec3(x / len, y / len, z / len);
			}
			return *this;
		}

		void normalize() {
			T len = length();
			if (len > 0) {
				x /= len;
				y /= len;
				z /= len;
			}
		}

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

		// Comparison operators
		constexpr bool operator==(const Vec3& v) const {
			return x == v.x && y == v.y && z == v.z;
		}

		constexpr bool operator!=(const Vec3& v) const {
			return !(*this == v);
		}
	};

	template<typename T>
	constexpr Vec2<T> operator*(T scalar, const Vec2<T>& v) {
		return v * scalar;
	}

	template<typename T>
	constexpr Vec3<T> operator*(T scalar, const Vec3<T>& v) {
		return v * scalar;
	}

	using Vec2f = Vec2<float>;
	using Vec2d = Vec2<double>;
	using Vec2i = Vec2<int>;

	using Vec3f = Vec3<float>;
	using Vec3d = Vec3<double>;
	using Vec3i = Vec3<int>;

	template<typename T>
	std::ostream& operator<<(std::ostream& os, const Vec2<T>& v) {
		return os << "Vec2(" << v.x << ", " << v.y << ")";
	}

	template<typename T>
	std::ostream& operator<<(std::ostream& os, const Vec3<T>& v) {
		return os << "Vec3(" << v.x << ", " << v.y << ", " << v.z << ")";
	}

	Vec3f barycentric(Vec2i* pts, Vec2i P);
}