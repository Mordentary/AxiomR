#pragma once
#include <cmath>
#include <iostream>
#include <type_traits>
#include <array>
#include <cassert>
#include <numbers>
#include <typeinfo>

#define GLM_ENABLE_EXPERIMENTAL
#define GLM_FORCE_INLINE
#define GLM_FORCE_SSE2
#define GLM_FORCE_AVX
#include "glm\glm.hpp"
#include <glm/gtx/fast_square_root.hpp>
namespace AR {
	//Color
	struct Color {
		unsigned char r, g, b, a;
	};
#define KB(x) ((x) * 1024ULL)
#define MB(x) ((x) * 1024ULL * 1024ULL)
	//	//============================================================
	//	// Vec2: 2D Vector
	//	//============================================================
	//	template<typename T>
	//	class Vec2 {
	//		static_assert(std::is_arithmetic<T>::value, "Vec2 requires an arithmetic type.");
	//	public:
	//		T x, y;
	//
	//		// Constructors
	//		constexpr Vec2() : x(T(0)), y(T(0)) {}
	//		constexpr Vec2(T xVal, T yVal) : x(xVal), y(yVal) {}
	//		constexpr Vec2(const Vec2& v) = default;
	//		Vec2& operator=(const Vec2& v) = default;
	//
	//		// Dot product
	//		constexpr T dot(const Vec2& v) const {
	//			return x * v.x + y * v.y;
	//		}
	//
	//		// Squared length
	//		constexpr T lengthSquared() const {
	//			return x * x + y * y;
	//		}
	//
	//		// Length (magnitude)
	//		T length() const {
	//			return std::sqrt(lengthSquared());
	//		}
	//
	//		// Normalized vector (returns a new vector)
	//		Vec2 normalized() const {
	//			T len = length();
	//			if (len > T(0)) {
	//				return Vec2(x / len, y / len);
	//			}
	//			// If length is zero, return the same vector
	//			return *this;
	//		}
	//
	//		// Normalize this vector in-place
	//		void normalize() {
	//			T len = length();
	//			if (len > T(0)) {
	//				x /= len;
	//				y /= len;
	//			}
	//		}
	//
	//		// Arithmetic operators
	//		constexpr Vec2 operator+(const Vec2& v) const {
	//			return Vec2(x + v.x, y + v.y);
	//		}
	//
	//		constexpr Vec2 operator-(const Vec2& v) const {
	//			return Vec2(x - v.x, y - v.y);
	//		}
	//
	//		constexpr Vec2 operator*(T scalar) const {
	//			return Vec2(x * scalar, y * scalar);
	//		}
	//
	//		constexpr Vec2 operator/(T scalar) const {
	//			return Vec2(x / scalar, y / scalar);
	//		}
	//
	//		constexpr Vec2& operator+=(const Vec2& v) {
	//			x += v.x;
	//			y += v.y;
	//			return *this;
	//		}
	//
	//		constexpr Vec2& operator-=(const Vec2& v) {
	//			x -= v.x;
	//			y -= v.y;
	//			return *this;
	//		}
	//
	//		constexpr Vec2& operator*=(T scalar) {
	//			x *= scalar;
	//			y *= scalar;
	//			return *this;
	//		}
	//
	//		constexpr Vec2& operator/=(T scalar) {
	//			x /= scalar;
	//			y /= scalar;
	//			return *this;
	//		}
	//
	//		constexpr bool operator==(const Vec2& v) const {
	//			return (x == v.x) && (y == v.y);
	//		}
	//
	//		constexpr bool operator!=(const Vec2& v) const {
	//			return !(*this == v);
	//		}
	//	};
	//
	//	// Scalar multiplication (scalar * vector)
	//	template<typename T>
	//	constexpr Vec2<T> operator*(T scalar, const Vec2<T>& v) {
	//		return v * scalar;
	//	}
	//
	//	//============================================================
	//	// Vec3: 3D Vector
	//	//============================================================
	//	template<typename T>
	//	class Vec3 {
	//		static_assert(std::is_arithmetic<T>::value, "Vec3 requires an arithmetic type.");
	//	public:
	//		T x, y, z;
	//
	//		// Constructors
	//		constexpr Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
	//		constexpr Vec3(T xVal, T yVal, T zVal) : x(xVal), y(yVal), z(zVal) {}
	//		constexpr Vec3(T val) : x(val), y(val), z(val) {}
	//		constexpr Vec3(const Vec3& v) = default;
	//		Vec3& operator=(const Vec3& v) = default;
	//
	//		// Dot product
	//		constexpr T dot(const Vec3& v) const {
	//			return x * v.x + y * v.y + z * v.z;
	//		}
	//
	//		constexpr Vec3 pow(const Vec3& v) const {
	//			return Vec3(
	//				std::pow(x, v.x),
	//				std::pow(y, v.y),
	//				std::pow(z, v.z)
	//			);
	//		}
	//		constexpr Vec3 cross(const Vec3& v) const {
	//			return Vec3(
	//				(y * v.z - z * v.y),
	//				(z * v.x - x * v.z),
	//				(x * v.y - y * v.x)
	//			);
	//		}
	//		constexpr Vec3 cross(Vec3&& v) const {
	//			return Vec3(
	//				(y * v.z - z * v.y),
	//				(z * v.x - x * v.z),
	//				(x * v.y - y * v.x)
	//			);
	//		}
	//
	//		// Squared length
	//		constexpr T lengthSquared() const {
	//			return x * x + y * y + z * z;
	//		}
	//
	//		// Length (magnitude)
	//		T length() const {
	//			return std::sqrt(lengthSquared());
	//		}
	//
	//		// Normalized vector (returns a new vector)
	//		Vec3 normalized() const {
	//			T len = length();
	//			if (len > T(0)) {
	//				return Vec3(x / len, y / len, z / len);
	//			}
	//			// If length is zero, return the same vector
	//			return *this;
	//		}
	//
	//		// Normalize this vector in-place
	//		void normalize() {
	//			T len = length();
	//			if (len > T(0)) {
	//				x /= len;
	//				y /= len;
	//				z /= len;
	//			}
	//		}
	//
	//		// Arithmetic operators
	//		constexpr Vec3 operator+(const Vec3& v) const {
	//			return Vec3(x + v.x, y + v.y, z + v.z);
	//		}
	//
	//		constexpr Vec3 operator-(const Vec3& v) const {
	//			return Vec3(x - v.x, y - v.y, z - v.z);
	//		}
	//
	//		constexpr Vec3 operator*(const Vec3& v) const {
	//			return Vec3(x * v.x, y * v.y, z * v.z);
	//		}
	//		constexpr Vec3 operator/(const Vec3& v) const {
	//			return Vec3(x / v.x, y / v.y, z / v.z);
	//		}
	//
	//		constexpr Vec3 operator-(T scalar) const {
	//			return Vec3(x - scalar, y - scalar, z - scalar);
	//		}
	//
	//		constexpr Vec3 operator*(T scalar) const {
	//			return Vec3(x * scalar, y * scalar, z * scalar);
	//		}
	//
	//		constexpr Vec3 operator/(T scalar) const {
	//			return Vec3(x / scalar, y / scalar, z / scalar);
	//		}
	//
	//		constexpr Vec3 operator-() const {
	//			return Vec3(-x, -y, -z);
	//		}
	//
	//		constexpr Vec3& operator+=(const Vec3& v) {
	//			x += v.x;
	//			y += v.y;
	//			z += v.z;
	//			return *this;
	//		}
	//
	//		constexpr Vec3& operator-=(const Vec3& v) {
	//			x -= v.x;
	//			y -= v.y;
	//			z -= v.z;
	//			return *this;
	//		}
	//
	//		constexpr Vec3& operator*=(T scalar) {
	//			x *= scalar;
	//			y *= scalar;
	//			z *= scalar;
	//			return *this;
	//		}
	//
	//		constexpr Vec3& operator/=(T scalar) {
	//			x /= scalar;
	//			y /= scalar;
	//			z /= scalar;
	//			return *this;
	//		}
	//
	//		constexpr bool operator==(const Vec3& v) const {
	//			return (x == v.x) && (y == v.y) && (z == v.z);
	//		}
	//
	//		constexpr bool operator!=(const Vec3& v) const {
	//			return !(*this == v);
	//		}
	//	};
	//
	//	// Scalar multiplication (scalar * vector)
	//	template<typename T>
	//	constexpr Vec3<T> operator*(T scalar, const Vec3<T>& v) {
	//		return v * scalar;
	//	}
	//
	//	//============================================================
	//	// Vec4: 4D Vector
	//	//============================================================
	//	template<typename T>
	//	class Vec4 {
	//		static_assert(std::is_arithmetic<T>::value, "Vec4 requires an arithmetic type.");
	//	public:
	//		T x, y, z, w;
	//
	//		// Constructors
	//		constexpr Vec4() : x(T(0)), y(T(0)), z(T(0)), w(T(0)) {}
	//		constexpr Vec4(T xVal, T yVal, T zVal, T wVal) : x(xVal), y(yVal), z(zVal), w(wVal) {}
	//		constexpr Vec4(T val) : x(val), y(val), z(val), w(val) {}
	//		constexpr Vec4(const Vec4& v) = default;
	//		Vec4& operator=(const Vec4& v) = default;
	//
	//		// Dot product
	//		constexpr T dot(const Vec4& v) const {
	//			return x * v.x + y * v.y + z * v.z + w * v.w;
	//		}
	//
	//		// Squared length
	//		constexpr T lengthSquared() const {
	//			return x * x + y * y + z * z + w * w;
	//		}
	//
	//		// Length (magnitude)
	//		T length() const {
	//			return std::sqrt(lengthSquared());
	//		}
	//
	//		// Normalized vector (returns a new vector)
	//		Vec4 normalized() const {
	//			T len = length();
	//			if (len > T(0)) {
	//				return Vec4(x / len, y / len, z / len, w / len);
	//			}
	//			// If length is zero, return the same vector
	//			return *this;
	//		}
	//
	//		// Normalize this vector in-place
	//		void normalize() {
	//			T len = length();
	//			if (len > T(0)) {
	//				x /= len;
	//				y /= len;
	//				z /= len;
	//				w /= len;
	//			}
	//		}
	//
	//		// Arithmetic operators
	//		constexpr Vec4 operator+(const Vec4& v) const {
	//			return Vec4(x + v.x, y + v.y, z + v.z, w + v.w);
	//		}
	//
	//		constexpr Vec4 operator-(const Vec4& v) const {
	//			return Vec4(x - v.x, y - v.y, z - v.z, w - v.w);
	//		}
	//
	//		constexpr Vec4 operator*(T scalar) const {
	//			return Vec4(x * scalar, y * scalar, z * scalar, w * scalar);
	//		}
	//
	//		constexpr Vec4 operator/(T scalar) const {
	//			return Vec4(x / scalar, y / scalar, z / scalar, w / scalar);
	//		}
	//
	//		constexpr Vec4 operator-() const {
	//			return Vec4(-x, -y, -z, -w);
	//		}
	//
	//		constexpr Vec4& operator+=(const Vec4& v) {
	//			x += v.x;
	//			y += v.y;
	//			z += v.z;
	//			w += v.w;
	//			return *this;
	//		}
	//
	//		constexpr Vec4& operator-=(const Vec4& v) {
	//			x -= v.x;
	//			y -= v.y;
	//			z -= v.z;
	//			w -= v.w;
	//			return *this;
	//		}
	//
	//		constexpr Vec4& operator*=(T scalar) {
	//			x *= scalar;
	//			y *= scalar;
	//			z *= scalar;
	//			w *= scalar;
	//			return *this;
	//		}
	//
	//		constexpr Vec4& operator/=(T scalar) {
	//			x /= scalar;
	//			y /= scalar;
	//			z /= scalar;
	//			w /= scalar;
	//			return *this;
	//		}
	//
	//		constexpr bool operator==(const Vec4& v) const {
	//			return (x == v.x) && (y == v.y) && (z == v.z) && (w == v.w);
	//		}
	//
	//		constexpr bool operator!=(const Vec4& v) const {
	//			return !(*this == v);
	//		}
	//	};
	//
	//	// Scalar multiplication (scalar * vector)
	//	template<typename T>
	//	constexpr Vec4<T> operator*(T scalar, const Vec4<T>& v) {
	//		return v * scalar;
	//	}
	//
	//	//============================================================
	//	// Type Aliases for Vectors
	//	//============================================================
	//	using glm::vec2 = Vec2<float>;
	//	using Vec2d = Vec2<double>;
	//	using glm::uvec2 = Vec2<int>;
	//	using Vec2u = Vec2<uint32_t>;
	//
	//	using glm::vec3 = Vec3<float>;
	//	using Vec3d = Vec3<double>;
	//	using Vec3i = Vec3<int>;
	//
	//	using glm::vec4 = Vec4<float>;
	//	using Vec4d = Vec4<double>;
	//	using Vec4i = Vec4<int>;
	//
	//	//============================================================
	//	// I/O Stream Operators for Vectors
	//	//============================================================
	//	template<typename T>
	//	std::ostream& operator<<(std::ostream& os, const Vec2<T>& v) {
	//		return os << "Vec2(" << v.x << ", " << v.y << ")";
	//	}
	//
	//	template<typename T>
	//	std::ostream& operator<<(std::ostream& os, const Vec3<T>& v) {
	//		return os << "Vec3(" << v.x << ", " << v.y << ", " << v.z << ")";
	//	}
	//
	//	template<typename T>
	//	std::ostream& operator<<(std::ostream& os, const Vec4<T>& v) {
	//		return os << "Vec4(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")";
	//	}
	//
	//	//============================================================
	//	   // Mat: Template Matrix Class
	//	   //============================================================
	//	template<size_t Rows, size_t Cols, typename T>
	//	class Mat {
	//		static_assert(Rows > 0 && Cols > 0, "Matrix must have positive dimensions.");
	//		static_assert(std::is_arithmetic<T>::value, "Matrix requires an arithmetic type.");
	//	public:
	//		std::array<std::array<T, Cols>, Rows> m;
	//
	//		// Constructors
	//		constexpr Mat() {
	//			for (size_t r = 0; r < Rows; ++r)
	//				for (size_t c = 0; c < Cols; ++c)
	//					m[r][c] = T(0);
	//		}
	//
	//		constexpr Mat(const std::array<std::array<T, Cols>, Rows>& values) : m(values) {}
	//
	//		// Access operators
	//		constexpr std::array<T, Cols>& operator[](size_t row) {
	//			assert(row < Rows);
	//			return m[row];
	//		}
	//
	//		constexpr const std::array<T, Cols>& operator[](size_t row) const {
	//			assert(row < Rows);
	//			return m[row];
	//		}
	//
	//		// Identity matrix (only for square matrices)
	//		template<size_t R = Rows, size_t C = Cols>
	//		static typename std::enable_if<R == C, Mat<R, C, T>>::type identity() {
	//			Mat<R, C, T> I;
	//			for (size_t r = 0; r < R; ++r)
	//				for (size_t c = 0; c < C; ++c)
	//					I.m[r][c] = (r == c) ? T(1) : T(0);
	//			return I;
	//		}
	//
	//		// Matrix multiplication
	//		template<size_t OtherCols>
	//		Mat<Rows, OtherCols, T> operator*(const Mat<Cols, OtherCols, T>& other) const {
	//			Mat<Rows, OtherCols, T> result;
	//			for (size_t r = 0; r < Rows; ++r) {
	//				for (size_t c = 0; c < OtherCols; ++c) {
	//					T sum = T(0);
	//					for (size_t k = 0; k < Cols; ++k) {
	//						sum += m[r][k] * other.m[k][c];
	//					}
	//					result.m[r][c] = sum;
	//				}
	//			}
	//			return result;
	//		}
	//		// Transpose of a 3x3 matrix
	//		Mat<3, 3, T> transpose() const {
	//			Mat<3, 3, T> result;
	//			for (size_t i = 0; i < 3; ++i) {
	//				for (size_t j = 0; j < 3; ++j) {
	//					result.m[i][j] = m[j][i];
	//				}
	//			}
	//			return result;
	//		}
	//
	//		// Inverse of a 3x3 matrix (you might want to add a check for singularity)
	//		Mat<3, 3, T> inverse() const {
	//			Mat<3, 3, T> result;
	//
	//			// Calculate the determinant of the 3x3 matrix
	//			T det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
	//				m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
	//				m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
	//
	//			// Adjugate matrix and then divide by determinant
	//			result.m[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det;
	//			result.m[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det;
	//			result.m[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det;
	//			result.m[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det;
	//			result.m[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det;
	//			result.m[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det;
	//			result.m[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det;
	//			result.m[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / det;
	//			result.m[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det;
	//
	//			return result;
	//		}
	//		// Matrix2x3-Vector3 multiplication
	//		template<typename U = T>
	//		typename std::enable_if<Rows == 2 && Cols == 3, Vec2<U>>::type
	//			operator*(const Vec3<U>& v) const {
	//			Vec2<U> result;
	//			result.x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z;
	//			result.y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z;
	//			return result;
	//		}
	//
	//		// Matrix-Vector multiplication (Vec2)
	//		template<typename U = T>
	//		typename std::enable_if<Rows == 2 && Cols == 2, Vec2<U>>::type
	//			operator*(const Vec2<U>& v) const {
	//			Vec2<U> result;
	//			result.x = m[0][0] * v.x + m[0][1] * v.y;
	//			result.y = m[1][0] * v.x + m[1][1] * v.y;
	//			return result;
	//		}
	//
	//		// Matrix-Vector multiplication (Vec3)
	//		template<typename U = T>
	//		typename std::enable_if<Rows == 3 && Cols == 3, Vec3<U>>::type
	//			operator*(const Vec3<U>& v) const {
	//			Vec3<U> result;
	//			result.x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z;
	//			result.y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z;
	//			result.z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z;
	//			return result;
	//		}
	//
	//		// Matrix-Vector multiplication (Vec4)
	//		template<typename U = T>
	//		typename std::enable_if<Rows == 4 && Cols == 4, Vec4<U>>::type
	//			operator*(const Vec4<U>& v) const {
	//			Vec4<U> result;
	//			result.x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3] * v.w;
	//			result.y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3] * v.w;
	//			result.z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3] * v.w;
	//			result.w = m[3][0] * v.x + m[3][1] * v.y + m[3][2] * v.z + m[3][3] * v.w;
	//			return result;
	//		}
	//
	//		// Addition
	//		Mat<Rows, Cols, T> operator+(const Mat<Rows, Cols, T>& other) const {
	//			Mat<Rows, Cols, T> result;
	//			for (size_t r = 0; r < Rows; ++r)
	//				for (size_t c = 0; c < Cols; ++c)
	//					result.m[r][c] = m[r][c] + other.m[r][c];
	//			return result;
	//		}
	//
	//		// Subtraction
	//		Mat<Rows, Cols, T> operator-(const Mat<Rows, Cols, T>& other) const {
	//			Mat<Rows, Cols, T> result;
	//			for (size_t r = 0; r < Rows; ++r)
	//				for (size_t c = 0; c < Cols; ++c)
	//					result.m[r][c] = m[r][c] - other.m[r][c];
	//			return result;
	//		}
	//
	//		// Scalar multiplication
	//		Mat<Rows, Cols, T> operator*(T scalar) const {
	//			Mat<Rows, Cols, T> result;
	//			for (size_t r = 0; r < Rows; ++r)
	//				for (size_t c = 0; c < Cols; ++c)
	//					result.m[r][c] = m[r][c] * scalar;
	//			return result;
	//		}
	//
	//		// Scalar division
	//		Mat<Rows, Cols, T> operator/(T scalar) const {
	//			Mat<Rows, Cols, T> result;
	//			for (size_t r = 0; r < Rows; ++r)
	//				for (size_t c = 0; c < Cols; ++c)
	//					result.m[r][c] = m[r][c] / scalar;
	//			return result;
	//		}
	//
	//		// In-place addition
	//		Mat<Rows, Cols, T>& operator+=(const Mat<Rows, Cols, T>& other) {
	//			for (size_t r = 0; r < Rows; ++r)
	//				for (size_t c = 0; c < Cols; ++c)
	//					m[r][c] += other.m[r][c];
	//			return *this;
	//		}
	//
	//		// In-place subtraction
	//		Mat<Rows, Cols, T>& operator-=(const Mat<Rows, Cols, T>& other) {
	//			for (size_t r = 0; r < Rows; ++r)
	//				for (size_t c = 0; c < Cols; ++c)
	//					m[r][c] -= other.m[r][c];
	//			return *this;
	//		}
	//
	//		// In-place scalar multiplication
	//		Mat<Rows, Cols, T>& operator*=(T scalar) {
	//			for (size_t r = 0; r < Rows; ++r)
	//				for (size_t c = 0; c < Cols; ++c)
	//					m[r][c] *= scalar;
	//			return *this;
	//		}
	//
	//		// In-place scalar division
	//		Mat<Rows, Cols, T>& operator/=(T scalar) {
	//			for (size_t r = 0; r < Rows; ++r)
	//				for (size_t c = 0; c < Cols; ++c)
	//					m[r][c] /= scalar;
	//			return *this;
	//		}
	//
	//		// Equality
	//		bool operator==(const Mat<Rows, Cols, T>& other) const {
	//			for (size_t r = 0; r < Rows; ++r)
	//				for (size_t c = 0; c < Cols; ++c)
	//					if (m[r][c] != other.m[r][c])
	//						return false;
	//			return true;
	//		}
	//
	//		// Inequality
	//		bool operator!=(const Mat<Rows, Cols, T>& other) const {
	//			return !(*this == other);
	//		}
	//
	//		//============================================================
	//		// set_col and get_col Functions
	//		//============================================================
	//
	//		// get_col for 2-row matrices
	//		template<size_t R = Rows>
	//		typename std::enable_if<R == 2, Vec2<T>>::type
	//			get_col(size_t col) const {
	//			assert(col < Cols && "Column index out of range.");
	//			return Vec2<T>(m[0][col], m[1][col]);
	//		}
	//
	//		// get_col for 3-row matrices
	//		template<size_t R = Rows>
	//		typename std::enable_if<R == 3, Vec3<T>>::type
	//			get_col(size_t col) const {
	//			assert(col < Cols && "Column index out of range.");
	//			return Vec3<T>(m[0][col], m[1][col], m[2][col]);
	//		}
	//
	//		// get_col for 4-row matrices
	//		template<size_t R = Rows>
	//		typename std::enable_if<R == 4, Vec4<T>>::type
	//			get_col(size_t col) const {
	//			assert(col < Cols && "Column index out of range.");
	//			return Vec4<T>(m[0][col], m[1][col], m[2][col], m[3][col]);
	//		}
	//
	//		// set_col for 2-row matrices
	//		template<size_t R = Rows>
	//		typename std::enable_if<R == 2, void>::type
	//			set_col(size_t col, const Vec2<T>& vec) {
	//			assert(col < Cols && "Column index out of range.");
	//			m[0][col] = vec.x;
	//			m[1][col] = vec.y;
	//		}
	//
	//		// set_col for 3-row matrices
	//		template<size_t R = Rows>
	//		typename std::enable_if<R == 3, void>::type
	//			set_col(size_t col, const Vec3<T>& vec) {
	//			assert(col < Cols && "Column index out of range.");
	//			m[0][col] = vec.x;
	//			m[1][col] = vec.y;
	//			m[2][col] = vec.z;
	//		}
	//
	//		// set_col for 4-row matrices
	//		template<size_t R = Rows>
	//		typename std::enable_if<R == 4, void>::type
	//			set_col(size_t col, const Vec4<T>& vec) {
	//			assert(col < Cols && "Column index out of range.");
	//			m[0][col] = vec.x;
	//			m[1][col] = vec.y;
	//			m[2][col] = vec.z;
	//			m[3][col] = vec.w;
	//		}
	//
	//		//============================================================
	//		// Static Methods for Specific Matrix Types (e.g., 4x4)
	//		//============================================================
	//
	//		// Translation Matrix (4x4)
	//		static Mat<4, 4, T> translate(const Vec3<T>& t) {
	//			static_assert(Rows == 4 && Cols == 4, "Translation is defined for 4x4 matrices.");
	//			Mat<4, 4, T> M = Mat<4, 4, T>::identity();
	//			M.m[0][3] = t.x;
	//			M.m[1][3] = t.y;
	//			M.m[2][3] = t.z;
	//			return M;
	//		}
	//
	//		// Rotation around X-axis (4x4)
	//		static Mat<4, 4, T> rotateX(T angleRadians) {
	//			static_assert(Rows == 4 && Cols == 4, "Rotation is defined for 4x4 matrices.");
	//			Mat<4, 4, T> R = Mat<4, 4, T>::identity();
	//			T c = std::cos(angleRadians);
	//			T s = std::sin(angleRadians);
	//			R.m[1][1] = c;
	//			R.m[1][2] = -s;
	//			R.m[2][1] = s;
	//			R.m[2][2] = c;
	//			return R;
	//		}
	//
	//		// Rotation around Y-axis (4x4)
	//		static Mat<4, 4, T> rotateY(T angleRadians) {
	//			static_assert(Rows == 4 && Cols == 4, "Rotation is defined for 4x4 matrices.");
	//			Mat<4, 4, T> R = Mat<4, 4, T>::identity();
	//			T c = std::cos(angleRadians);
	//			T s = std::sin(angleRadians);
	//			R.m[0][0] = c;
	//			R.m[0][2] = s;
	//			R.m[2][0] = -s;
	//			R.m[2][2] = c;
	//			return R;
	//		}
	//
	//		// Rotation around Z-axis (4x4)
	//		static Mat<4, 4, T> rotateZ(T angleRadians) {
	//			static_assert(Rows == 4 && Cols == 4, "Rotation is defined for 4x4 matrices.");
	//			Mat<4, 4, T> R = Mat<4, 4, T>::identity();
	//			T c = std::cos(angleRadians);
	//			T s = std::sin(angleRadians);
	//			R.m[0][0] = c;
	//			R.m[0][1] = -s;
	//			R.m[1][0] = s;
	//			R.m[1][1] = c;
	//			return R;
	//		}
	//
	//		// Combined Rotation (4x4)
	//		static Mat<4, 4, T> rotateXYZ(T angleX, T angleY, T angleZ) {
	//			Mat<4, 4, T> rx = rotateX(angleX);
	//			Mat<4, 4, T> ry = rotateY(angleY);
	//			Mat<4, 4, T> rz = rotateZ(angleZ);
	//			return rz * ry * rx;
	//		}
	//
	//		// LookAt Matrix (4x4)
	//		static Mat<4, 4, T> lookAt(const Vec3<T>& eye, const Vec3<T>& target, const Vec3<T>& up) {
	//			static_assert(Rows == 4 && Cols == 4, "LookAt is defined for 4x4 matrices.");
	//			Vec3<T> f = (target - eye).normalized();
	//			Vec3<T> r = f.cross(up).normalized();
	//			Vec3<T> u = r.cross(f);
	//
	//			Mat<4, 4, T> view = Mat<4, 4, T>::identity();
	//			view.m[0][0] = r.x; view.m[0][1] = r.y; view.m[0][2] = r.z; view.m[0][3] = -r.dot(eye);
	//			view.m[1][0] = u.x; view.m[1][1] = u.y; view.m[1][2] = u.z; view.m[1][3] = -u.dot(eye);
	//			view.m[2][0] = -f.x; view.m[2][1] = -f.y; view.m[2][2] = -f.z; view.m[2][3] = f.dot(eye);
	//			return view;
	//		}
	//
	//		// Perspective Projection Matrix (4x4)
	//		static Mat<4, 4, T> perspective(T fovRadians, T aspect, T nearZ, T farZ) {
	//			static_assert(Rows == 4 && Cols == 4, "Perspective is defined for 4x4 matrices.");
	//			Mat<4, 4, T> P;
	//			T f = T(1) / std::tan(fovRadians / T(2));
	//			P.m[0][0] = f / aspect;
	//			P.m[1][1] = f;
	//			P.m[2][2] = (farZ + nearZ) / (nearZ - farZ);
	//			P.m[2][3] = (T(2) * farZ * nearZ) / (nearZ - farZ);
	//			P.m[3][2] = -T(1);
	//			P.m[3][3] = T(0);
	//			return P;
	//		}
	//	};
	//
	//	//============================================================
	//	// Type Aliases for Common Matrices
	//	//============================================================
	//	template<typename T>
	//	using mat2 = Mat<2, 2, T>;
	//
	//	template<typename T>
	//	using mat3 = Mat<3, 3, T>;
	//
	//	template<typename T>
	//	using mat4 = Mat<4, 4, T>;
	//
	//	using mat2f = mat2<float>;
	//	using mat2d = mat2<double>;
	//	using mat2i = mat2<int>;
	//
	//	using glm::mat3 = mat3<float>;
	//	using mat3d = mat3<double>;
	//	using mat3i = mat3<int>;
	//
	//	using glm::mat4 = mat4<float>;
	//	using mat4d = mat4<double>;
	//	using mat4i = mat4<int>;
	//
	//	//============================================================
	//	// I/O Stream Operators for Matrices
	//	//============================================================
	//	template<size_t Rows, size_t Cols, typename T>
	//	std::ostream& operator<<(std::ostream& os, const Mat<Rows, Cols, T>& mat) {
	//		os << "Mat<" << Rows << ", " << Cols << ", " << typeid(T).name() << ">{\n";
	//		for (size_t r = 0; r < Rows; ++r) {
	//			os << "  ";
	//			for (size_t c = 0; c < Cols; ++c) {
	//				os << mat.m[r][c];
	//				if (c < Cols - 1) os << ", ";
	//			}
	//			os << "\n";
	//		}
	//		os << "}";
	//		return os;
	//	}
	//
	//	//============================================================
	//	// Additional Utility Functions
	//	//============================================================
	//	inline glm::vec4 toVec4f(const glm::vec3& vec3, float w = 1.0f) {
	//		return glm::vec4(vec3.x, vec3.y, vec3.z, w);
	//	}
	//	inline glm::vec3 perspectiveDivision(const glm::vec4& vec4) {
	//		if (vec4.w != 0.0f) {
	//			return glm::vec3(vec4.x / vec4.w, vec4.y / vec4.w, vec4.z / vec4.w);
	//		}
	//		// If w is zero, return the vector as is (could represent a direction)
	//		return glm::vec3(vec4.x, vec4.y, vec4.z);
	//	}
	//
	//	inline glm::vec3 toVec3f(const glm::vec4& vec4) {
	//		if (vec4.w != 0.0f) {
	//			return glm::vec3(vec4.x, vec4.y, vec4.z);
	//		}
	//		// If w is zero, return the vector as is (could represent a direction)
	//		return glm::vec3(vec4.x, vec4.y, vec4.z);
	//	}
	//	double degreeToRad(double degrees);
	//	// Barycentric Coordinates (Assuming 2D Points)
	//	glm::vec3 barycentric(const std::array<glm::uvec2, 3>& pts, const glm::uvec2& P);
	//
	//	glm::vec3 reflect(const glm::vec3& incident, const glm::vec3& normal);
}