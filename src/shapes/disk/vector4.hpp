#ifndef __VECTOR4__
#define __VECTOR4__

#include <cmath>
#include <cstdint>
#include <iostream>

namespace v4l {

	template<class T>
	class vector4
	{
	public:
		constexpr vector4(T n = T(0)) : x(n), y(n), z(n), o(n) {}
		constexpr vector4(T x, T y, T z, T o) : x(x), y(y), z(z), o(o) {}

		constexpr vector4(vector4<T>& v) : x(v.x), y(v.y), z(v.z), o(v.o) {}
		constexpr vector4(const vector4<T>& v) : x(v.x), y(v.y), z(v.z), o(v.o) {}

		~vector4() {}
	public:
		T x, y, z, o;
	};

	// Defines

	using vector4d = vector4<double>;
	using vector4f = vector4<float>;

	// Scalar operations

	template<class T> constexpr inline
		vector4<T> operator*(const vector4<T>& v, const T& n)
	{
		return vector4<T>(v.x * n, v.y * n, v.z * n, v.o * n);
	}

	template<class T> constexpr inline
		vector4<T> operator*(const T& n, const vector4<T>& v)
	{
		return vector4<T>(n * v.x, n * v.y, n * v.z, n * v.o);
	}

	template<class T> constexpr inline
		vector4<T> operator/(const vector4<T>& v, const T& n)
	{
		return vector4<T>(v.x / n, v.y / n, v.z / n, v.o / n);
	}

	template<class T> constexpr inline
		vector4<T> operator/(const T& n, const vector4<T>& v)
	{
		return vector4<T>(n / v.x, n / v.y, n / v.z, n / v.o);
	}

	// Element-wise operations

	template<class T> constexpr inline
		vector4<T> operator-(const vector4<T>& v)
	{
		return vector4<T>(-v.x, -v.y, -v.z, -v.o);
	}

	template<class T> constexpr inline
		vector4<T> operator+(const vector4<T>& v1, const vector4<T>& v2)
	{
		return vector4<T>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z, v1.o + v2.o);
	}

	template<class T> constexpr inline
		vector4<T> operator-(const vector4<T>& v1, const vector4<T>& v2)
	{
		return vector4<T>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.o - v2.o);
	}

	template<class T> constexpr inline
		vector4<T> operator*(const vector4<T>& v1, const vector4<T>& v2)
	{
		return vector4<T>(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z, v1.o * v2.o);
	}

	template<class T> constexpr inline
		vector4<T> operator/(const vector4<T>& v1, const vector4<T>& v2)
	{
		return vector4<T>(v1.x / v2.x, v1.y / v2.y, v1.z / v2.z, v1.o / v2.o);
	}

	// Vector operations

	template<class T> constexpr inline
		T dot(const vector4<T>& v1, const vector4<T>& v2)
	{
		return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.o * v2.o;
	}

	template<class T> constexpr inline
		T operator|(const vector4<T>& v1, const vector4<T>& v2)
	{
		return dot(v1, v2);
	}

	template<class T> constexpr inline
		T length(const vector4<T>& v)
	{
		return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z + v.o * v.o);
	}

	template<class T> constexpr inline
		T length2(const vector4<T>& v)
	{
		return v.x * v.x + v.y * v.y + v.z * v.z + v.o * v.o;
	}

	template<class T> constexpr inline
		vector4<T> normalize(const vector4<T>& v)
	{
		T norm = T(1) / length(v);
		return vector4<T>(v.x * norm, v.y * norm, v.z * norm, v.o * norm);
	}

	template<class T> inline
		std::ostream & operator<<(std::ostream & os, const vector4<T>& v)
	{
		os << "(" << v.x << ", " << v.y << ", " << v.z << ", " << v.o << ")";
		return os;
	}

} // namespace v3l

#endif // __VECTOR3__
