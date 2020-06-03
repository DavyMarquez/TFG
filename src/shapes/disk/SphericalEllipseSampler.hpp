/*
 * SphericalEllipseSampler.hpp
 *
 *      Author: ibon
 */

#ifndef __SPHERICALELLIPSESAMPLING_HPP__
#define __SPHERICALELLIPSESAMPLING_HPP__

#include <algorithm>
#include <numeric>
#include <iostream>

#define BOOST_MATH_MAX_SERIES_ITERATION_POLICY 500
#include <boost/math/special_functions/ellint_3.hpp>
namespace {
	using boost::math::ellint_3;
}
#include <boost/math/special_functions/ellint_1.hpp>
namespace {
	using boost::math::ellint_1;
}

#include "vector3.hpp"
namespace {
	using v3l::vector3;
}

#include "mathUtil.hpp"

// Select appropiate precomputed table size
//#include "CDFTable32x32.hpp"
//#include "CDFTable64x64.hpp"
//#include "CDFTable128x128.hpp"
//#include "CDFTable256x256.hpp"
#include "CDFTable512x512.hpp"
//#include "CDFTable1024x1024.hpp"


template<class T>
class SphericalEllipseSampler
{
public:
	enum class Type {
		ZERO, // Indicates unsampleable case
		POLAR,
		POLARTAB,
		CONCENTRIC,
		CONCENTRICTAB,
		CYLINDRICAL,
		CYLINDRICALTAB
	} type;

	// Parameters of each algorithm variant
	union {
		struct BoothData {
			// Elliptical integral constants
			T c, n, k;
			// Precalculated constants
			T phi_coeff, t_coeff, asqrtb, bsqrta;
			// Quadrant area
			T S;
		} boothData;
		struct {
			// Tabulated CDF
			double* CDF1D;
			// Quadrant area
			T S;
		} boothTabData;
		struct UrenaData {
			// Spherical ellipse arcs
			T alpha, beta;
			// Tangent ellipse params
			T a_p, b_p;
			// Elliptical integral constants
			T n, k, m;
			// Precalculated constants
			T p, c1;
			// Half area
			T S_2;
		} urenaData;
	} data;

	// Spherical ellipse params
	T a, b;
private:
	// Constants
	static constexpr T pi = T(3.14159265358979323846);
	static constexpr T pi_2 = T(3.14159265358979323846 / 2.0);
	static constexpr T pi_4 = T(3.14159265358979323846 / 4.0);
	static constexpr T two_pi = T(2.0*3.14159265358979323846);
	static constexpr T eps = T(1e-6);
	static constexpr T tol = T(1e-5);
public:
	vector3<T> sample(const T e1, const T e2, T& pdf) const;

	T samplePdf(const vector3<T>& p) const;
private:
	/*
	 * Azimuth angle that covers the specified fractional area over the spherical ellipse,
	 * as found by inverting Booth's expression by Newton iterations
	 * 		u   : fraction of total area
	 * 		pdf : probability of the generated sample
	 * 		tol : tolerance
	 */
	T BoothInversionNewton(const T u, T &pdf, const T tol = tol) const;
	/*
	 * Azimuth angle that covers the specified fractional area over the spherical ellipse,
	 * as found by inverting a tabulated approximation by spherical triangles
	 * 		u   : fraction of total area
	 * 		pdf : probability of the generated sample
	 * 		h_i : height bound of the sampled spherical triangle
	 */
	T BoothInversionTable(const T u, T& pdf, T& h_i) const;
	/*
	 * Probability of sampling a given azimuth angle using table inversion
	 * 		phi : sampled azimuth angle
	 */
	T BoothPdfTable(const T phi) const;
	/*
	 * Fractional minor arc that covers the specified fractional area over the spherical ellipse,
	 * as found by inverting Urena's expression by Newton iterations
	 * 		u   : fraction of total area
	 * 		tol : tolerance
	 * 		pdf : probability of the generated sample
	 */
	T UrenaInversionNewton(T u, T &pdf, T tol = tol) const;
};

template<class T> inline
vector3<T> SphericalEllipseSampler<T>::sample(const T e1, const T e2, T& pdf) const
{
	vector3<T> p(T(0));

	switch (this->type) {
	case Type::POLAR:
	case Type::POLARTAB:
	{
		// This is a continuous variation of polar mapping that has a single discontinuity
		// and avoids the disk-to-square-and-back step. The branching could probably be
		// greatly reduced

		// Choose sampling quadrant
		T u;
		if (e1 >= T(0.75)) {
			// Fourth quadrant
			u = T(1) - T(4) * (e1 - T(0.75));
		}
		else if (e1 >= T(0.5)) {
			// Third quadrant
			u = T(4) * (e1 - T(0.5));
		}
		else if (e1 >= T(0.25)) {
			// Second quadrant
			u = T(1) - T(4) * (e1 - T(0.25));
		}
		else {
			// First quadrant
			u = T(4) * e1;
		}

		T v = e2;

		// Sample azimuth angle
		T phi_u, h_i = T(0);
		if (this->type == Type::POLAR) {
			// Newton root finding
			phi_u = this->BoothInversionNewton(u, pdf);
		}
		else {
			// Tabulated CDF inversion
			phi_u = this->BoothInversionTable(u, pdf, h_i);
		}

		// Translate to sampling quadrant
		if (e1 >= T(0.75)) {
			// Fourth quadrant
			phi_u = two_pi - phi_u;
		}
		else if (e1 >= T(0.5)) {
			// Third quadrant
			phi_u = pi + phi_u;
		}
		else if (e1 >= T(0.25)) {
			// Second quadrant
			phi_u = pi - phi_u;
		}

		const T a = this->a;
		const T b = this->b;

		// Calculate spherical radius and height for the point over the
		// spherical ellipse edge defined by phi
		T sinphi = std::sin(phi_u);
		T cosphi = std::cos(phi_u);

		T r_u = (a*b) / std::sqrt(a*a * sinphi*sinphi + b * b * cosphi*cosphi);
		T h_u = std::sqrt(T(1) - r_u * r_u);

		T h_vu;
		if (this->type == Type::POLAR) {
			// Sample uniformly a point along the spherical ellipse arc
			h_vu = T(1) - (T(1) - h_u) * v;
		}
		else {
			// Sample uniformly a point along the spherical triangle arc
			h_vu = T(1) - (T(1) - h_i) * v;

			// Check if the sampled point is inside the spherical ellipse
			// or reject otherwise
			if (h_vu < h_u) {
				pdf = T(0);
				return vector3<T>(T(0));
			}
		}

		T r_vu = std::sqrt(1 - h_vu * h_vu);

		// Local coordinates
		T x_vu = r_vu * cosphi;
		T y_vu = r_vu * sinphi;
		T z_vu = h_vu;

		p = vector3<T>(x_vu, y_vu, z_vu);
	}
	break;
	case Type::CONCENTRIC:
	case Type::CONCENTRICTAB:
	{
		// A variant of Shirley's clever continuous mapping that unfortunately contains tons
		// of branching because we can only sample the first quadrant of the spherical
		// ellipse, so we need to go back and forth between quadrants

		T u = T(2) * e1 - T(1);
		T v = T(2) * e2 - T(1);

		// Choose sampling quadrant
		T r, theta;
		if (u == T(0) && v == T(0)) {
			r = T(0);
			theta = T(0);
		}
		else if (std::abs(u) > std::abs(v)) {
			r = u;
			if (u < T(0)) {
				// Third quadrant
				theta = pi + pi_4 * v / u;
			}
			else if (v >= T(0)) {
				// First quadrant positive
				theta = pi_4 * v / u;
			}
			else {
				// First quadrant negative
				theta = two_pi + pi_4 * v / u;
			}
		}
		else {
			r = v;
			if (v >= T(0)) {
				// Second quadrant
				theta = pi_2 - pi_4 * u / v;
			}
			else {
				// Fourth quadrant
				theta = 3 * pi_2 - pi_4 * u / v;
			}
		}

		// Translate sampling coords
		T u_p;
		if (theta > 3 * pi_2) {
			u_p = 1 - (theta - 3 * pi_2) / pi_2;
		}
		else if (theta > pi) {
			u_p = (theta - pi) / pi_2;
		}
		else if (theta > pi_2) {
			u_p = 1 - (theta - pi_2) / pi_2;
		}
		else {
			u_p = theta / pi_2;
		}

		T v_p = r * r;

		// Sample azimuth angle
		T phi, h_i = 0.0;
		if (this->type == Type::CONCENTRIC) {
			// Newton root finding
			phi = this->BoothInversionNewton(u_p, pdf);
		}
		else {
			// Tabulated CDF inversion
			phi = this->BoothInversionTable(u_p, pdf, h_i);
		}

		// Translate sampled angle to sampling quadrant
		T phi_u;
		if (theta > 3 * pi_2) {
			phi_u = two_pi - phi;
		}
		else if (theta > pi) {
			phi_u = pi + phi;
		}
		else if (theta > pi_2) {
			phi_u = pi - phi;
		}
		else {
			phi_u = phi;
		}

		const T a = this->a;
		const T b = this->b;

		// Calculate spherical radius and height for the point over the spherical ellipse edge
		// delimited by angle phi
		T sinphi = std::sin(phi_u);
		T cosphi = std::cos(phi_u);

		T r_u = (a*b) / std::sqrt(a*a * sinphi*sinphi + b * b * cosphi*cosphi);
		T h_u = std::sqrt(T(1) - r_u * r_u);

		T h_vu;
		if (this->type == Type::CONCENTRIC) {
			// Sample uniformly a point along the spherical ellipse arc
			h_vu = T(1) - (T(1) - h_u) * v_p;
		}
		else {
			// Sample uniformly a point along the spherical triangle arc
			h_vu = T(1) - (T(1) - h_i) * v_p;

			// Check if the sampled point is inside the spherical ellipse
			// or reject otherwise
			if (h_vu < h_u) {
				pdf = T(-1);
				return vector3<T>(T(0));
			}
		}

		T r_vu = std::sqrt(T(1) - h_vu * h_vu);

		// Local coordinates
		T x_vu = r_vu * cosphi;
		T y_vu = r_vu * sinphi;
		T z_vu = h_vu;

		p = vector3<T>(x_vu, y_vu, z_vu);
	}
	break;
	case Type::CYLINDRICAL:
	case Type::CYLINDRICALTAB:
	{
		T beta_t = this->UrenaInversionNewton(e1, pdf);

		const T a_p = this->data.urenaData.a_p;
		const T b_p = this->data.urenaData.b_p;

		// Aboslute coords. of the two points that delimit the tangent
		// ellipse's chord defined by beta
		T y_1 = std::tan(beta_t);
		T x_1 = a_p * std::sqrt(T(1) - (y_1 * y_1) / (b_p * b_p));

		// Limiting points reprojected into the sphere
//			vector3<T>s0_p = normalize(vector3<T>(-x_1, y_1 , T(1)));
		vector3<T>s1_p = normalize(vector3<T>(x_1, y_1, T(1)));

		// Linear interpolation between limiting points
		T w = T(2) * e2 - T(1);
		T f = std::sqrt((T(1) - std::pow(w * s1_p.x, 2)) /
			(T(1) - std::pow(s1_p.x, 2)));

		p = vector3<T>(s1_p.x * w, s1_p.y * f, s1_p.z * f);
	}
	break;
	case Type::ZERO:
	{
		pdf = T(0);
		p = vector3<T>(T(0));
	}
	}

	return p;
}

//TODO: Find pdf by Table Inversion
template<class T> inline
T SphericalEllipseSampler<T>::samplePdf(const vector3<T>&) const
{
	switch (this->type) {
	case Type::POLAR:
	case Type::CONCENTRIC:
	{
		auto& data = this->data.boothData;
		return T(1) / (T(4) * data.S);
	}
	break;
	case Type::POLARTAB:
	case Type::CONCENTRICTAB:
	{
		auto& data = this->data.boothTabData;
		return T(1) / (T(4) * data.S);
	}
	break;
	case Type::CYLINDRICAL:
	case Type::CYLINDRICALTAB:
	{
		auto& data = this->data.urenaData;
		return T(1) / (T(2) * data.S_2);
	}
	break;
	case Type::ZERO:
	{
		return T(0);
	}
	}

	return T(0);
}

template<class T> inline
T SphericalEllipseSampler<T>::BoothInversionNewton(const T u, T &pdf, const T tol) const
{
	auto& data = this->data.boothData;

	// Avoid going outside of function domain
	T phi_max = u * T(pi_2);

	// Initial guess, assume planar ellipse
	// TODO: Better initial approximation
	T phi = std::atan(data.phi_coeff * std::tan(phi_max));
	phi = sml::clamp(phi, T(0), phi_max);

	// Target fractional area
	T S_u = u * data.S;

	// Newton iterations
	for (size_t i = 0; i < 100; i++) {
		T sinphi = std::sin(phi);

		T sinphi2 = sinphi * sinphi;
		T cosphi2 = T(1) - sinphi2;

		// Parametric angle
		T t = sml::clamp(std::atan(data.t_coeff * std::tan(phi)), T(0), T(pi_2));

		T sint = std::sin(t);
		T sint2 = sint * sint;

		// Surface area for phi
		T omega = phi - data.c * ellint_3(data.k, data.n, t);

		T diff = omega - S_u;

		// Check if required tolerance has been reached
		if (std::abs(diff) < tol) {
			break;
		}

		// Derivative at phi
		T der = T(1) - ((data.c * data.asqrtb * data.bsqrta) /
			((T(1) - data.n * sint2) * std::sqrt(T(1) - data.k*data.k * sint2) *
			((data.bsqrta) * cosphi2 + (data.asqrtb) * sinphi2)));

		// Avoid nearly zero divisions
		if (std::abs(der) < T(eps)) {
			break;
		}

		// Avoid going outside of function domain
		phi = sml::clamp(phi - diff / der, T(0), phi_max);
	}

	pdf = T(1) / (T(4) * data.S);

	return phi;
}

template<class T> inline
T SphericalEllipseSampler<T>::BoothInversionTable(const T u, T& pdf, T& h_i) const
{
	auto& data = this->data.boothTabData;

	// Binary search
	double* ptr = std::upper_bound(data.CDF1D, data.CDF1D + CDF_SIZE, u);
	size_t index = size_t(std::max(ptrdiff_t(0), ptr - data.CDF1D));

	T CDF_0 = (index > 0) ? T(data.CDF1D[index - 1]) : T(0);
	T CDF_1 = T(data.CDF1D[index]);
	T CDF_D = CDF_1 - CDF_0;

	// Calculate height at the start of the table entry
	T phi_i = (T(index) * T(pi_2)) / T(CDF_SIZE);
	T sinphi_i = std::sin(phi_i);
	T cosphi_i = std::cos(phi_i);

	const T a = this->a;
	const T b = this->b;

	T rphi_i = (a*b) / std::sqrt(a*a * sinphi_i*sinphi_i + b * b * cosphi_i*cosphi_i);
	h_i = std::sqrt(T(1) - rphi_i * rphi_i);

	// Find interpolated azimuth angle for sampling
	T du = (u - CDF_0) / CDF_D;
	T phi = ((T(index) + du) * T(pi_2)) / T(CDF_SIZE);

	pdf = T(1) / (T(4) * data.S);

	return phi;
}

template<class T> inline
T SphericalEllipseSampler<T>::UrenaInversionNewton(T u, T &pdf, T tol) const
{
	auto& data = this->data.urenaData;

	// Initial guess, assume linear relation
	T t = u;

	// Target fractional area
	T S_u = T(2) * data.S_2 * u;

	// Newton iterations
	for (size_t i = 0; i < 100; i++) {
		T beta_t = std::abs(data.beta * (T(2) * t - T(1)));

		T Omega_t;
		if (t == T(0)) {
			Omega_t = T(0);
		}
		else if (t == T(0.5)) {
			Omega_t = data.S_2;
		}
		else if (t == T(1)) {
			Omega_t = T(2) * data.S_2;
		}
		else {
			T ampl = -std::asin(std::min(T(1), std::tan(beta_t) / data.b_p));
			T I = data.c1 * data.b_p * (data.p * ellint_1(data.k, ampl) -
				(data.p + T(1)) * ellint_3(data.k, data.n, ampl));

			// Partial area
			Omega_t = data.S_2 + sml::sign(t - T(0.5)) * I;
		}

		T diff = Omega_t - S_u;

		// Check if required tolerance is reached
		if (std::abs(diff) < tol) {
			break;
		}

		T tan_beta_t = std::tan(beta_t);
		T delta = data.c1 * std::sqrt((T(1) - data.p * tan_beta_t * tan_beta_t) /
			(T(1) - data.m * data.p * tan_beta_t * tan_beta_t));

		T Omega_t_d = T(2) * data.beta * delta;

		// Avoid nearly zero divisions
		if (std::abs(Omega_t_d) < T(eps)) {
			break;
		}

		t = t - (Omega_t - S_u) / Omega_t_d;
	}

	pdf = T(1) / (T(2) * data.S_2);

	return (T(2) * t - T(1)) * data.beta;
}

#endif // __SPHERICALELLIPSESAMPLING_HPP__
