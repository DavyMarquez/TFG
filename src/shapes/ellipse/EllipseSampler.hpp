#ifndef __ELLIPSESAMPLER_HPP__
#define __ELLIPSESAMPLER_HPP__

#include <iostream>
#include <fstream>
#include <mitsuba/core/matrix.h>
MTS_NAMESPACE_BEGIN

#define BOOST_MATH_MAX_SERIES_ITERATION_POLICY 500
//#include <boost/math/special_functions/ellint_3.hpp>
#include "../dependencies/include/boost/math/special_functions/ellint_3.hpp"


namespace {
	using boost::math::ellint_3;
}
//#include <boost/math/special_functions/ellint_1.hpp>
#include "../dependencies/include/boost/math/special_functions/ellint_1.hpp"
namespace {
	using boost::math::ellint_1;
}
#include "vector3.hpp"
namespace {
	using v3l::vector3;
}

#include "SphericalEllipseSampler.hpp"

template<class T>
class EllipseSampler;

template<class T>
class EllipseSamplingRecord
{
	friend class EllipseSampler<T>;
private:
	// Type of sampling
	SphericalEllipseSampler<T> algorithm;
	// Ellipse coordinate system
	vector3<T> o, xd, yd, zd;
	// Distance to ellipse plane from shading point
	T d;
	// Reprojection constants
	T yhp_norm, dhp_norm;
public:
	EllipseSamplingRecord() {};
	/*
	 * A point over the unit sphere corresponding to the mapping from the [0..1]^2
	 * square to the ellipse projected solid angle
	 * 		u   : x coordinate over the [0..1]^2 square
	 * 		v   : y coordinate over the [0..1]^2 square
	 * 		pdf : probability of the sample
	 */
	vector3<T> sample(const T e1, const T e2, T& pdf) const;
	/*
	 * pdf corresponding to the ellipse sample as generated by this method
	 */
	T samplePdf(const vector3<T>& p) const;
};

template<class T> inline
vector3<T> EllipseSamplingRecord<T>::sample(const T e1, const T e2, T& pdf) const
{
	// Generate sample in local coords of the spherical ellipse
	vector3<T> local = this->algorithm.sample(e1, e2, pdf);

	// Don't reproject discarded samples
	if (pdf == T(-1)) {
		return vector3<T>(T(0));
	}

	T xu = local.x, yu = local.y, zu = local.z;

	// Ellipse's local coords over the sphere
	T yu_p = yu * this->dhp_norm + zu * this->yhp_norm;
	T zu_p = -yu * this->yhp_norm + zu * this->dhp_norm;

	T prop = this->d / zu_p;

	// Ellipse's local coords projected on the ellipse
	T x = xu * prop;
	T y = yu_p * prop;
	T z = this->d;

	// Global coords over the ellipse
	/*vector3<T> result = this->o + x * this->xd + y * this->yd + z * this->zd
	cout << "global coords over the ellipse" << endl;
	cout << result.x << ", " << result.y << ", " << result.z << endl;*/
	return this->o + x * this->xd + y * this->yd + z * this->zd;
}

template<class T> inline
T EllipseSamplingRecord<T>::samplePdf(const vector3<T>& p) const
{
	return this->algorithm.samplePdf(p);
}

// -----------------------------------------------------------------------------

template<class T>
class EllipseSampler
{
private:
	// Type of sampling
	typename SphericalEllipseSampler<T>::Type algorithm;

	// Tabulated CDF
	double* CDF;
private:
	// Constants
	static constexpr T pi = T(3.14159265358979323846);
	static constexpr T pi_2 = T(3.14159265358979323846 / 2.0);
	static constexpr T eps = T(1e-20);
public:
	EllipseSampler() {}

	EllipseSampler(const typename SphericalEllipseSampler<T>::Type algo) :
		algorithm(algo), CDF(nullptr)
	{
		// Load CDF data
		this->CDF = (double *)&CDF_hex_data;
	}

	/*
	 * Ellipse sampling record
	 * 		o  : shading point
	 * 		semiaxisA  : ellipse major semiaxis (scaled)
	 *		semiaxisB  : ellipse minor semiaxis (scaled)
	 *      xd, yd, zd, cd : ellipse coordinate system
	 *      
	 */
	EllipseSamplingRecord<T> createRecord(vector3<T>& o, 
					const T semiaxisA, const T semiaxisB, 
					const vector3<T>& xd,
					const vector3<T>& yd,
					const vector3<T>& zd,
					const vector3<T>& cd) const;
};

template<class T> inline
EllipseSamplingRecord<T> EllipseSampler<T>::createRecord(vector3<T>& o, 
						const T semiaxisA, 
						const T semiaxisB, 
						const vector3<T>& xd,
						const vector3<T>& yd,
						const vector3<T>& zd,
						const vector3<T>& cd) const
{

	EllipseSamplingRecord<T> record;
	/*for (int px = 100; px >= 0; px--) {
		for (int py = 100; py >= 0; py--) {
			for(int pz = 100; pz >= 0; pz--){*/
				// Arbitrary ellipse semiaxes nearly zero
				if (semiaxisA < eps || semiaxisB < eps) {
					record.algorithm.type = SphericalEllipseSampler<T>::Type::ZERO;
					return record;
				}

				/*o.x = px;
				o.y = py;
				o.z = pz;*/
				//record.o = v3l::vector3<T>(0, 4.9, 0);
				//shading point to arbitrary ellipse center
				vector3<T> l = cd - o;

				//Distance from shading point to arbitrary ellipse plane
				T d = -v3l::dot(l, zd);

				// Check if the arbitrary ellipse normal is pointing to the opposite side
				// and calculate ellipse's coordinates Z axis
				if (d < T(0)) {
					record.zd = zd / length(zd);
					d = -d;
				}
				else {
					record.zd = -zd / length(zd);
				}
				record.d = d;

				//Point inside arbitrary ellipse plane
				if (d == 0) {
					//cout << "punto dentro del plano de la elipse" << endl;
					record.algorithm.type = SphericalEllipseSampler<T>::Type::ZERO;
					return record;
				}

				// Build reference system
				if (d < T(0)) {
					record.xd = -xd / length(xd);
					record.yd = -yd / length(yd);
				}
				else {
					record.xd = xd / length(xd);
					record.yd = yd / length(yd);
				}

				//Store local coordinate system origin
				record.o = o;
				//record.o = v3l::vector3<T>(px, py, pz);
				//record.o = v3l::vector3<T>(0,0,5.5);
				/* GET THE TANGENT ELLIPSE USING THE HEITZ TRANSFORM METHOD
				 * https://hal.archives-ouvertes.fr/hal-01561624/document */

				TVector3<double> shadPoint = TVector3<double>(o.x, o.y, o.z); //Shading point
				//TVector3<double> shadPoint = TVector3<double>(-1.0, -1.0, -1.0);
				//TVector3<double> shadPoint = TVector3<double>(px, py, pz);
				//TVector3<double> shadPoint = TVector3<double>(0, 0, 5.5);
				//cout << "shadPoint = " << shadPoint.toString() << endl;
				
				TVector3<double> C = TVector3<double>(cd.x, cd.y, cd.z); //Ellipse center
				T lx = semiaxisA, //Both ellipse semiaxes
					ly = semiaxisB;

				//Shading point transform matrix
				double toWorldMShadAux[4][4] = { {1.0, 0.0, 0.0, shadPoint.x},
											{0.0, 1.0, 0.0, shadPoint.y},
											{0.0, 0.0, 1.0, shadPoint.z},
											{0.0, 0.0, 0.0, 1.0}};

				Matrix<4, 4, double> toWorldMShad = Matrix<4,4,double>(toWorldMShadAux);
				Matrix<4,4,double> toWorldMShad_inv;
				bool useless = toWorldMShad.invert(toWorldMShad_inv);

				//Arbitrary ellipse transform matrix
				double toWorldMrealAux[4][4] = { {xd.x, yd.x, zd.x, cd.x},
									{xd.y, yd.y, zd.y, cd.y},
									{xd.z, yd.z, zd.z, cd.z},
									{0.0, 0.0, 0.0, 1.0} };
				Matrix<4, 4, double> toWorldMreal = Matrix<4,4,double>(toWorldMrealAux);

				//Arbitrary ellipse transform matrix over shading point
				Matrix<4, 4, double> toWorldM = toWorldMShad_inv * toWorldMreal;
				C = TVector3<double>(toWorldM(0, 3), toWorldM(1, 3), toWorldM(2, 3));
				//cout << std::setprecision(5);
				
				//Arbitrary ellipse coordinate system with center in {0,0,0}
				double matrixAux[4][4] = {
					{toWorldM(0,0),toWorldM(0,1),toWorldM(0,2),0},
					{toWorldM(1,0),toWorldM(1,1),toWorldM(1,2),0},
					{toWorldM(2,0),toWorldM(2,1),toWorldM(2,2),0},
					{toWorldM(3,0),toWorldM(3,1),toWorldM(3,2),0} };
				
				Matrix<4, 4, double> matrix = Matrix<4, 4, double>(matrixAux);

				//Arbitrary ellipse axes
				TVector3<double> Vx = TVector3<double>(toWorldM(0,0), toWorldM(1,0), toWorldM(2,0)),
					Vy = TVector3<double>(toWorldM(0, 1), toWorldM(1, 1), toWorldM(2, 1)),
					Vz = TVector3<double>(toWorldM(0, 2), toWorldM(1, 2), toWorldM(2, 2));
				Vx = Vx / Vx.length();
				Vy = Vy / Vy.length();
				Vz = Vz / Vz.length();

				//Arbitrary ellipse center in this basis
				double xc = C.x * Vx.x + C.y * Vx.y + C.z * Vx.z,
					yc = C.x * Vy.x + C.y * Vy.y + C.z * Vy.z,
					zc = C.x * Vz.x + C.y * Vz.y + C.z * Vz.z;

				double Qaux[3][3] = { {(zc*zc) / (lx*lx), 0, (-zc / (lx*lx))*xc},
				{0, (zc*zc) / (ly*ly), (-zc / (ly*ly))*yc},
				{(-zc / (lx*lx))*xc, (-zc / (ly*ly))*yc,
					((xc*xc) / (lx*lx)) + ((yc*yc) / (ly*ly)) - 1} };

				Matrix<3, 3, double> Q = Matrix<3, 3, double>(Qaux);
				
				/*cout << "Q2 = " << Q2 << endl;
				cout << "Q = " << Q.toString() << endl;*/

				//eigendecomposition
				double E[3]; //eigenvalues
				Q.symEig(Q, E); //Q now stores the eigenvectors

				Matrix<3, 3, double> Q2 = Matrix<3,3, double>(cos(4));

				/*double bo = Q2(0, 0);
				cout << std::setprecision(15) << bo << endl;*/
				/*char varr;
				std::cin >> varr;*/

				int order[3] = { 0, 1, 2 };
				//order Q by descending eigenvalues
				for (int i = 0; i < 3 - 1; i++) {
					int maxIndex = i;
					for (int j = i + 1; j < 3; j++) {
						if (E[j] > E[maxIndex]) {
							maxIndex = j;
						}
					}
					double aux = E[maxIndex];
					E[maxIndex] = E[i];
					E[i] = aux;
					int auxI = order[maxIndex];
					order[maxIndex] = order[i];
					order[i] = auxI;
				}
				double Qaux2[3][3] = { {Q(0,order[0]), Q(0,order[1]), Q(0,order[2])},
				{Q(1,order[0]), Q(1,order[1]), Q(1,order[2])},
				{Q(2,order[0]), Q(2,order[1]), Q(2,order[2])} };

				Q = Matrix<3,3,double>(Qaux2);

				/*scaling factor so that tangent ellipse center distance to shading point
				* equals 1 */
				double lambda = 1;

				//Tangent ellipse axes
				double semiAxisA_tangent = lambda * sqrt(-E[2] / E[0]),
					semiAxisB_tangent = lambda * sqrt(-E[2] / E[1]);

				TVector3<double> center = TVector3<double>(lambda * Q(0,2), lambda * Q(1, 2), lambda * Q(2, 2)),
					center2 = TVector3<double>(-lambda * Q(0, 2), -lambda * Q(1, 2), -lambda * Q(2, 2));

				// "_h" refers to homogenoeus coordinates
				double centerMatrixAux[4][1] = { center[0], center[1], center[2], 1.0 };
				double centerMatrixAux2[4][1] = { center2[0], center2[1], center2[2], 1.0 };
				Matrix<4, 1, double> centerAux_h = matrix * Matrix<4,1, double>(centerMatrixAux),
					center2Aux_h = matrix * Matrix<4, 1, double>(centerMatrixAux2);

				TVector3<double> centerAux = TVector3<double>(centerAux_h(0,0), centerAux_h(0,1), centerAux_h(0,2)),
					center2Aux = TVector3<double>(center2Aux_h(0,0), center2Aux_h(0,1), center2Aux_h(0,2));

				//keep tangent ellipse center closer to arbitrary ellipse
				if ((centerAux - C).length() > (center2Aux - C).length()) {
					center = center2;
					centerAux = center2Aux;
					centerAux_h = center2Aux_h;
				}

				TVector3<double> Xaxis = TVector3<double>(Q(0,0), Q(1, 0), Q(2, 0)),
					Yaxis = TVector3<double>(Q(0, 1), Q(1, 1), Q(2, 1)),
					Zaxis = TVector3<double>(Q(0, 2), Q(1, 2), Q(2, 2));
				Xaxis = Xaxis / Xaxis.length();
				Yaxis = Yaxis / Yaxis.length();
				Zaxis = Zaxis / Zaxis.length();

				//Tangent ellipse global coordinate system 
				double Xaxis_hAux_[4][1] = { Xaxis.x, Xaxis.y, Xaxis.z, 0.0 };
				Matrix<4, 1, double> Xaxis_hAux = Matrix<4, 1, double>(Xaxis_hAux_);
				Matrix<4, 1, double> col1 = matrix * Xaxis_hAux;
				//col1 = matrix * Vector4(Xaxis.x, Xaxis.y, Xaxis.z, 0);
				//col2 = matrix * Vector4(Yaxis.x, Yaxis.y, Yaxis.z, 0);
				double Yaxis_hAux_[4][1] = { Yaxis.x, Yaxis.y, Yaxis.z, 0.0 };
				Matrix<4, 1, double> Yaxis_hAux = Matrix<4, 1, double>(Yaxis_hAux_);
				Matrix<4, 1, double> col2 = matrix * Yaxis_hAux;
				//col3 = matrix * Vector4(Zaxis.x, Zaxis.y, Zaxis.z, 0);
				double Zaxis_hAux_[4][1] = { Zaxis.x, Zaxis.y, Zaxis.z, 0.0 };
				Matrix<4, 1, double> Zaxis_hAux = Matrix<4, 1, double>(Zaxis_hAux_);
				Matrix<4, 1, double> col3 = matrix * Zaxis_hAux;
				//col4 = Vector4(centerAux.x, centerAux.y, centerAux.z, 1);
				double centerAux_hAux_[4][1] = { centerAux.x, centerAux.y, centerAux.z, 1.0 };
				Matrix<4, 1, double> centerAux_hAux = Matrix<4, 1, double>(centerAux_hAux_);
				Matrix<4, 1, double> col4 = centerAux_hAux;
				double toWorldM2Aux[4][4] = { {col1(0,0), col2(0,0), col3(0,0), col4(0,0)},
				{col1(0,1), col2(0,1), col3(0,1), col4(0,1)},
				{col1(0,2), col2(0,2), col3(0,2), col4(0,2)},
				{col1(0,3), col2(0,3), col3(0,3), col4(0,3)} };
				Matrix<4,4,double> toWorldM2 = Matrix<4,4,double>(toWorldM2Aux);
				toWorldM2 = toWorldMShad * toWorldM2;

				//Tangent ellipse axes in global coordinates
				Xaxis = TVector3<double>(toWorldM2(0, 0), toWorldM2(1,0), toWorldM2(2,0));
				Yaxis = TVector3<double>(toWorldM2(0, 1), toWorldM2(1, 2), toWorldM2(2, 1));
				Zaxis = TVector3<double>(toWorldM2(0, 2), toWorldM2(1, 2), toWorldM2(2, 2));

				//tangent ellipse center in global coordinates
				centerAux = TVector3<double>(toWorldM2(0,3), toWorldM2(1, 3), toWorldM2(2, 3));

				/*cout << "distance to tangent ellipse" << endl;
				cout << (centerAux - shadPoint).length() << endl;*/
				//Tangent ellipse semiaxes
				double at = semiAxisA_tangent,
					bt = semiAxisB_tangent;

				double otroAuxiliar[4][1] = { toWorldM2(0,3), toWorldM2(1, 3), toWorldM2(2, 3), 1.0 };
				Matrix<4, 1, double> centerAux_hmatrix = Matrix<4, 1, double>(otroAuxiliar);
				Matrix<4,1,double> centerAux_inv = toWorldMShad_inv * centerAux_hmatrix;
				//Spherical ellipse center coordinates
				/*double yh_p = centerAux_inv(0,0),
					dh_p = centerAux_inv(0,1);
				if (yh_p < double(0)) {
					yh_p = -yh_p;
				}
				if (dh_p < double(0)) {
					dh_p = -dh_p;
				}*/

				//cout << centerAux_inv.toString() << endl;
				/*cout << "shad point " << endl;
				cout << shadPoint.toString() << endl;*/
				//cout << "yh_p: " << yh_p << ", dh_p: " << dh_p << endl;
				//char variable;
				//std::cin >> variable;
				// Calculations in double precision to avoid numerical instabilities
				double alpha = atan(at),
					beta = atan(bt);

				//Spherical ellipse semiaxes
				double a = sin(alpha),
					b = sin(beta);
				if (a < b) {
					a = sin(beta);
					b = sin(alpha);
				}
				/*double a = at/ sqrt((at*at)+1), //problems with ellint_3 function
					b = bt / sqrt((bt*bt) + 1);
				if (a < b) {
					a = sqrt((bt*bt) + 1);
					b = sqrt((at*at) + 1);
				}*/

				//cout << "a: " << a << ", b: " << b << endl;
				/*cout << "alpha = " << alpha << ", beta = " << beta << endl;
				cout << "at = " << at << ", bt = " << bt << endl;*/
				// Discard very small ellipses
				// Possibly unnecessary when working with the tabulated variant

				/*cout << std::setprecision(15);
				cout << "shadPoint = " << shadPoint.toString() << endl;
				cout << "ELLIPSE: a = " << a << ", b = " << b << endl;*/

				/*char var;
				std::cin >> var;*/
				//cout << "a = " << a << ", b = " << b << endl;
				if (a < eps || b < eps) {
					/*cout << "semiaxes too small" << endl;
					char var;
					std::cin >> var;*/
					record.algorithm.type = SphericalEllipseSampler<T>::Type::ZERO;
					return record;
				}

				SphericalEllipseSampler<T> algo;
				switch (this->algorithm) {
				case SphericalEllipseSampler<T>::Type::POLAR:
				case SphericalEllipseSampler<T>::Type::CONCENTRIC:
				{
					//cout << "a = " << a << ", b = " << b << endl;
					double c = (b*(1.0 - a * a)) / (a * std::sqrt(1.0 - b * b));
					double k = std::sqrt((a*a - b * b) / (1.0 - b * b));
					double n = (a*a - b * b) / (a*a*(1.0 - b * b));
					//cout << "c: " << c << ", k: " << k << ", n: " << n << endl;
					double S = pi_2 - c * ellint_3(k, n, pi_2);
					T asqrtb = a * std::sqrt(T(1) - b * b);
					T bsqrta = b * std::sqrt(T(1) - a * a);
					T phi_coeff = (b*std::sqrt(T(1) - a * a)) / (a*std::sqrt(T(1) - b * b));
					T t_coeff = 1 / phi_coeff;
					/*cout << "S: " << S << ", asqrtb: " << asqrtb << ", bsqrta: " << bsqrta << endl;
					cout << "phi_coeff: " << phi_coeff << ", T_coeff: " << t_coeff << endl;
					char var;
					std::cin >> var;*/
					algo.type = this->algorithm;
					// Spherical ellipse axis
					algo.a = a;
					algo.b = b;
					// Elliptic integral constants
					algo.data.boothData.c = T(c);
					algo.data.boothData.n = T(n);
					algo.data.boothData.k = T(k);
					// Precalc. constants
					algo.data.boothData.phi_coeff = T(phi_coeff);
					algo.data.boothData.t_coeff = T(t_coeff);
					algo.data.boothData.asqrtb = T(asqrtb);
					algo.data.boothData.bsqrta = T(bsqrta);
					// Precalc quadrant area
					algo.data.boothData.S = S;
				}
				break;
				case SphericalEllipseSampler<T>::Type::POLARTAB:
				case SphericalEllipseSampler<T>::Type::CONCENTRICTAB:
				{
					// Find row of the precalc. CDF corresponding to the current arc proportion
					// TODO: Tabulate by (b/a) intervals instead of (beta/alpha) ones
					double alpha = std::asin(a), beta = std::asin(b);

					size_t index = size_t(std::floor(CDF_SIZE * (beta / alpha))) * CDF_SIZE;
					double* CDF1D = &this->CDF[index];

					double c = (b*(1.0 - a * a)) / (a * std::sqrt(1.0 - b * b));
					double k = std::sqrt((a*a - b * b) / (1.0 - b * b));
					double n = (a*a - b * b) / (a*a*(1.0 - b * b));

					double S = pi_2 - c * ellint_3(k, n, pi_2);

					algo.type = this->algorithm;
					// Spherical ellipse axis
					algo.a = a;
					algo.b = b;
					// 1D CDF for fast table inversion
					algo.data.boothTabData.CDF1D = CDF1D;
					// Precalc quadrant area
					algo.data.boothTabData.S = S;
				}
				break;
				case SphericalEllipseSampler<T>::Type::CYLINDRICAL:
				{
					T a_p = a / std::sqrt(T(1) - a * a); //a_p = tan(alpha)
					T b_p = b / std::sqrt(T(1) - b * b);

					T alpha = std::asin(a);
					T beta = std::asin(b);

					T p = T(1) / (b_p * b_p);
					T m = (a_p * a_p - b_p * b_p) / (a_p * a_p + 1);
					T k = std::sqrt(m);
					T n = -b_p * b_p;
					T c1 = (T(2) * a_p) / std::sqrt(T(1) + a_p * a_p);

					T S_2 = c1 * b_p * (p * ellint_1(k, -pi_2) - (p + T(1)) * ellint_3(k, n, -pi_2));

					algo.type = this->algorithm;
					// Spherical ellipse axis
					algo.a = a;
					algo.b = b;
					// Spherical ellipse arcs
					algo.data.urenaData.alpha = alpha;
					algo.data.urenaData.beta = beta;
					// Tangent ellipse params
					algo.data.urenaData.a_p = a_p;
					algo.data.urenaData.b_p = b_p;
					// Elliptic integral constants
					algo.data.urenaData.n = n;
					algo.data.urenaData.k = k;
					algo.data.urenaData.m = m;
					// Precalc. constants
					algo.data.urenaData.p = p;
					algo.data.urenaData.c1 = c1;
					// Precalc. half area
					algo.data.urenaData.S_2 = S_2;
				}
				break;
				case SphericalEllipseSampler<T>::Type::CYLINDRICALTAB:
				{
					// TODO: Some kind of unbiased tabulation scheme for cilindric sampling
				}
				break;
				case SphericalEllipseSampler<T>::Type::ZERO:
					break;
				}

				record.algorithm = algo;

				// Reprojection constants
				//T axis_norm = T(1) / std::sqrt(yh_p*yh_p + dh_p * dh_p);
				//record.yhp_norm = yh_p * axis_norm;

				//record.dhp_norm = dh_p * axis_norm;
				return record;
				/*cout << "shadPoint: {" << px << ", " << py << ", " << pz << "}" << endl;
				cout << "o: " << record.o << endl;
				cout << "xd: " << record.xd << endl;
				cout << "yd: " << record.yd << endl;
				cout << "zd: " << record.zd << endl;
				cout << "d: " << record.d << endl;
				cout << "yhp_norm: " << record.yhp_norm << endl;
				cout << "dhp_norm: " << record.dhp_norm << endl << endl;*/
			/*}
		}
	}
	char var;
	std::cin >> var;*/
	return record;
}

MTS_NAMESPACE_END
#endif // __ELLIPSESAMPLER_HPP__
