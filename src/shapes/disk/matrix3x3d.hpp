#ifndef __MATRIX3X3D__
#define __MATRIX3X3D__

#include <cmath>
#include <cstdint>
#include <iostream>

#include "vector3.hpp"
namespace {
	using v3l::vector3;
}

namespace m3d {

class matrix3x3d{
public:
	matrix3x3d(){
		m[0][0] = 0; m[0][1] = 0; m[0][2] = 0;
		m[1][0] = 0; m[1][1] = 0; m[1][2] = 0;
		m[2][0] = 0; m[2][1] = 0; m[2][2] = 0;
	} 
	matrix3x3d(double m00, double m01, double m02, 
					     double m10, double m11, double m12, 
					     double m20, double m21, double m22){
		m[0][0] = m00; m[0][1] = m01; m[0][2] = m02;
		m[1][0] = m10; m[1][1] = m11; m[1][2] = m12;
		m[2][0] = m20; m[2][1] = m21; m[2][2] = m22;
	} 

	matrix3x3d(v3l::vector3<double> col1, v3l::vector3<double> col2, 
			   v3l::vector3<double> col3){
		m[0][0] = col1.x; m[0][1] = col2.x; m[0][2] = col3.x;
		m[1][0] = col1.y; m[1][1] = col2.y; m[1][2] = col3.y;
		m[2][0] = col1.z; m[2][1] = col2.z; m[2][2] = col3.z;
	}

	v3l::vector3<double> col(int i){
		v3l::vector3<double> result = v3l::vector3<double>(m[0][i], m[1][i],
											m[2][i]);
		return result;
	}

	v3l::vector3<double> row(int i){
		v3l::vector3<double> result = v3l::vector3<double>(m[i][0], m[i][1],
											m[i][2]);
		return result;
	}

    void tred2(double d[3], double e[3]) {
        int N = 3;
        for (int j = 0; j < N; j++)
            d[j] = m[N - 1][j];

        // Householder reduction to tridiagonal form.
        for (int i = N - 1; i > 0; i--) {
            // Scale to avoid under/overflow.
            double scale = 0.0, h = 0.0;

            for (int k = 0; k < i; k++)
                scale = scale + std::abs(d[k]);

            if (scale == 0.0) {
                e[i] = d[i - 1];
                for (int j = 0; j < i; j++) {
                    d[j] = m[i - 1][j];
                    m[i][j] = 0.0;
                    m[j][i] = 0.0;
                }
            } else {
                // Generate Householder vector.

                for (int k = 0; k < i; k++) {
                    d[k] /= scale;
                    h += d[k] * d[k];
                }
                double f = d[i - 1],
                      g = std::sqrt(h);

                if (f > 0)
                    g = -g;
                e[i] = scale * g;
                h = h - f * g;
                d[i - 1] = f - g;

                for (int j = 0; j < i; j++)
                    e[j] = 0.0;

                 // Apply similarity transformation to remaining columns.
                for (int j = 0; j < i; j++) {
                    f = d[j];
                    m[j][i] = f;
                    g = e[j] + m[j][j] * f;
                    for (int k = j + 1; k <= i - 1; k++) {
                        g += m[k][j] * d[k];
                        e[k] += m[k][j] * f;
                    }
                    e[j] = g;
                }

                f = 0.0;
                for (int j = 0; j < i; j++) {
                    e[j] /= h;
                    f += e[j] * d[j];
                }
                double hh = f / (h + h);

                for (int j = 0; j < i; j++)
                    e[j] -= hh * d[j];
                for (int j = 0; j < i; j++) {
                    f = d[j];
                    g = e[j];
                    for (int k = j; k <= i - 1; k++)
                        m[k][j] -= (f * e[k] + g * d[k]);
                    d[j] = m[i - 1][j];
                    m[i][j] = 0.0;
                }
            }
            d[i] = h;
        }
    }

    inline void tql2(double d[3], double e[3]) {
        int N = 3;
        for (int i = 1; i < N; i++){
            e[i - 1] = e[i];
        }

        e[N - 1] = 0.0;
        double f = 0.0, tst1 = 0.0;

        double eps = pow(2.0, -52.0);

        for (int l = 0; l < N; l++) {
            // Find small subdiagonal element
            tst1 = std::max(tst1, std::abs(d[l]) + std::abs(e[l]));
            int mAux = l;

            while (mAux < N) {
                if (std::abs(e[mAux]) <= eps * tst1)
                    break;
                mAux++;
            }

            // If mAux == l, d[l] is an eigenvalue,
            // otherwise, iterate.
            if (mAux > l) {
                int iter = 0;

                do {
                    iter = iter + 1;    // (Could check iteration count here.)

                    // Compute implicit shift
                    double g = d[l];
                    double p = (d[l + 1] - g) / (2.0f * e[l]);
                    double r = math::hypot2((double) 1, p);

                    if (p < 0)
                        r = -r;
                    d[l] = e[l] / (p + r);
                    d[l + 1] = e[l] * (p + r);
                    double dl1 = d[l + 1];

                    double h = g - d[l];

                    for (int i = l + 2; i < N; i++)
                        d[i] -= h;
                    f = f + h;

                    // Implicit QL transformation.
                    p = d[mAux];
                    double c = 1.0f;
                    double c2 = c, c3 = c;
                    double el1 = e[l + 1];
                    double s = 0.0, s2 = 0.0;

                    for (int i = mAux - 1; i >= l; i--) {
                        c3 = c2;
                        c2 = c;
                        s2 = s;
                        g = c * e[i];
                        h = c * p;
                        r = math::hypot2(p, e[i]);
                        e[i + 1] = s * r;
                        s = e[i] / r;
                        c = p / r;
                        p = c * d[i] - s * g;
                        d[i + 1] = h + s * (c * g + s * d[i]);

                        // Accumulate transformation.
                        for (int k = 0; k < N; k++) {
                            h = m[k][i + 1];
                            m[k][i + 1] =
                                s * m[k][i] + c * h;
                            m[k][i] = c * m[k][i] - s * h;
                        }
                    }
                    p = -s * s2 * c3 * el1 * e[l] / dl1;
                    e[l] = s * p;
                    d[l] = c * p;
                    // Check for convergence.
                } while (std::abs(e[l]) > eps * tst1);
            }
            d[l] = d[l] + f;
            e[l] = 0.0;
        }

        // Sort eigenvalues and corresponding vectors.
        for (int i = 0; i < N - 1; i++) {
            int k = i;

            double p = d[i];

            for (int j = i + 1; j < N; j++) {
                if (d[j] < p) {
                    k = j;
                    p = d[j];
                }
            }

            if (k != i) {
                d[k] = d[i];
                d[i] = p;
                for (int j = 0; j < N; j++) {
                    p = m[j][i];
                    m[j][i] = m[j][k];
                    m[j][k] = p;
                }
            }
        }
    }

    inline void symEig(matrix3x3d& Q, double d[3]){
        double e[3];
        tred2(Q, e);
        tql2(d,e);
    }

	~matrix3x3d() {}
public:
	double m[3][3];
};


inline
std::ostream & operator<<(std::ostream & os, const matrix3x3d& v)
{
   os << "[[" << v.m[0][0] << ", " << v.m[0][1] << ", " << v.m[0][2] << "]" << endl;
   os << "[" << v.m[1][0] << ", " << v.m[1][1] << ", " << v.m[1][2] << "]" << endl;
   os << "[" << v.m[2][0] << ", " << v.m[2][1] << ", " << v.m[2][2] << "]" << endl;
   return os ;
}

inline
vector3<double> operator*(const matrix3x3d& m, const vector3<double>& v)
{
	return vector3<double>(m.m[0][0]*v.x + m.m[0][1]*v.y + m.m[0][2]*v.z,
		m.m[1][0]*v.x + m.m[1][1]*v.y + m.m[1][2]*v.z,
		m.m[2][0]*v.x + m.m[2][1]*v.y + m.m[2][2]*v.z);
}

inline
matrix3x3d operator*(const matrix3x3d& m1, const matrix3x3d& m2)
{
	matrix3x3d res = matrix3x3d();
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			res.m[i][j] = 0;
			for(int k = 0; k < 3; k++){
				res.m[i][j] += m1.m[i][k] * m2.m[k][j];
			}
		}
	}
	return res;
}

} // namespace m3d

#endif // __MAdoubleRIX3X3D__
