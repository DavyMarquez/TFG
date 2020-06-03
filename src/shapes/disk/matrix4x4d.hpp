/*
 * Beware, for its adapted to be used on coordinate system matrices
 */

#ifndef __MATRIX4X4D__
#define __MATRIX4X4D__

#include <cmath>
#include <cstdint>
#include <iostream>

#include "vector4.hpp"
namespace {
	using v4l::vector4;
}

namespace m4d {

class matrix4x4d{
public:
	matrix4x4d(){
		m[0][0] = 0; m[0][1] = 0; m[0][2] = 0; m[0][3] = 0;
		m[1][0] = 0; m[1][1] = 0; m[1][2] = 0; m[1][3] = 0;
		m[2][0] = 0; m[2][1] = 0; m[2][2] = 0; m[2][3] = 0;
		m[3][0] = 0; m[3][1] = 0; m[3][2] = 0; m[3][3] = 0;
	} 
	matrix4x4d(double m00, double m01, double m02, double m03, 
					     double m10, double m11, double m12, double m13, 
					     double m20, double m21, double m22, double m23, 
					     double m30, double m31, double m32, double m33){
		m[0][0] = m00; m[0][1] = m01; m[0][2] = m02; m[0][3] = m03;
		m[1][0] = m10; m[1][1] = m11; m[1][2] = m12; m[1][3] = m13;
		m[2][0] = m20; m[2][1] = m21; m[2][2] = m22; m[2][3] = m23;
		m[3][0] = m30; m[3][1] = m31; m[3][2] = m32; m[3][3] = m33;
	} 

	matrix4x4d(v4l::vector4<double> col1, v4l::vector4<double> col2, 
			   v4l::vector4<double> col3, v4l::vector4<double> col4){
		m[0][0] = col1.x; m[0][1] = col2.x; m[0][2] = col3.x; m[0][3] = col4.x;
		m[1][0] = col1.y; m[1][1] = col2.y; m[1][2] = col3.y; m[1][3] = col4.y;
		m[2][0] = col1.z; m[2][1] = col2.z; m[2][2] = col3.z; m[2][3] = col4.z;
		m[3][0] = col1.o; m[3][1] = col2.o; m[3][2] = col3.o; m[3][3] = col4.o;
	}

	v4l::vector4<double> col(int i){
		v4l::vector4<double> result = v4l::vector4<double>(m[0][i], m[1][i],
											m[2][i], m[3][i]);
		return result;
	}

	v4l::vector4<double> row(int i){
		v4l::vector4<double> result = v4l::vector4<double>(m[i][0], m[i][1],
											m[i][2], m[i][3]);
		return result;
	}

	~matrix4x4d() {}
public:
	double m[4][4];
};


inline
std::ostream & operator<<(std::ostream & os, const matrix4x4d& v)
{
   os << "[[" << v.m[0][0] << ", " << v.m[0][1] << ", " << v.m[0][2] << ", " << v.m[0][3] << "]" << endl;
   os << "[" << v.m[1][0] << ", " << v.m[1][1] << ", " << v.m[1][2] << ", " << v.m[1][3] << "]" << endl;
   os << "[" << v.m[2][0] << ", " << v.m[2][1] << ", " << v.m[2][2] << ", " << v.m[2][3] << "]" << endl;
   os << "[" << v.m[3][0] << ", " << v.m[3][1] << ", " << v.m[3][2] << ", " << v.m[3][3] << "]]" << endl;
   return os ;
}

inline
vector4<double> operator*(const matrix4x4d& m, const vector4<double>& v)
{
	return vector4<double>(m.m[0][0]*v.x + m.m[0][1]*v.y + m.m[0][2]*v.z + m.m[0][3]*v.o,
		m.m[1][0]*v.x + m.m[1][1]*v.y + m.m[1][2]*v.z + m.m[1][3]*v.o,
		m.m[2][0]*v.x + m.m[2][1]*v.y + m.m[2][2]*v.z + m.m[2][3]*v.o,
		m.m[3][0]*v.x + m.m[3][1]*v.y + m.m[3][2]*v.z + m.m[3][3]*v.o);
}

inline
matrix4x4d operator*(const matrix4x4d& m1, const matrix4x4d& m2)
{
	matrix4x4d res = matrix4x4d();
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			res.m[i][j] = 0;
			for(int k = 0; k < 4; k++){
				res.m[i][j] += m1.m[i][k] * m2.m[k][j];
			}
		}
	}
	return res;
}

inline
matrix4x4d invert(matrix4x4d m){
    matrix4x4d inv = matrix4x4d();
    double det;

    inv.m[0][0] = m.m[1][1]  * m.m[2][2] * m.m[3][3] - 
             m.m[1][1]  * m.m[2][3] * m.m[3][2] - 
             m.m[2][1]  * m.m[1][2]  * m.m[3][3] + 
             m.m[2][1]  * m.m[1][3]  * m.m[3][2] +
             m.m[3][1] * m.m[1][2]  * m.m[2][3] - 
             m.m[3][1] * m.m[1][3]  * m.m[2][2];

    inv.m[1][0] = -m.m[1][0]  * m.m[2][2] * m.m[3][3] + 
              m.m[1][0]  * m.m[2][3] * m.m[3][2] + 
              m.m[2][0]  * m.m[1][2]  * m.m[3][3] - 
              m.m[2][0]  * m.m[1][3]  * m.m[3][2] - 
              m.m[3][0] * m.m[1][2]  * m.m[2][3] + 
              m.m[3][0] * m.m[1][3]  * m.m[2][2];

    inv.m[2][0] = m.m[1][0]  * m.m[2][1] * m.m[3][3] - 
             m.m[1][0]  * m.m[2][3] * m.m[3][1] - 
             m.m[2][0]  * m.m[1][1] * m.m[3][3] + 
             m.m[2][0]  * m.m[1][3] * m.m[3][1] + 
             m.m[3][0] * m.m[1][1] * m.m[2][3] - 
             m.m[3][0] * m.m[1][3] * m.m[2][1];

    inv.m[3][0] = 0;
    /*inv.m[3][0] = -m.m[1][0]  * m.m[2][1] * m.m[3][2] + 
               m.m[1][0]  * m.m[2][2] * m.m[3][1] +
               m.m[2][0]  * m.m[1][1] * m.m[3][2] - 
               m.m[2][0]  * m.m[1][2] * m.m[3][1] - 
               m.m[3][0] * m.m[1][1] * m.m[2][2] + 
               m.m[3][0] * m.m[1][2] * m.m[2][1];*/

    inv.m[0][1] = -m.m[0][1]  * m.m[2][2] * m.m[3][3] + 
              m.m[0][1]  * m.m[2][3] * m.m[3][2] + 
              m.m[2][1]  * m.m[0][2] * m.m[3][3] - 
              m.m[2][1]  * m.m[0][3] * m.m[3][2] - 
              m.m[3][1] * m.m[0][2] * m.m[2][3] + 
              m.m[3][1] * m.m[0][3] * m.m[2][2];

    inv.m[1][1] = m.m[0][0]  * m.m[2][2] * m.m[3][3] - 
             m.m[0][0]  * m.m[2][3] * m.m[3][2] - 
             m.m[2][0]  * m.m[0][2] * m.m[3][3] + 
             m.m[2][0]  * m.m[0][3] * m.m[3][2] + 
             m.m[3][0] * m.m[0][2] * m.m[2][3] - 
             m.m[3][0] * m.m[0][3] * m.m[2][2];

    inv.m[2][1] = -m.m[0][0]  * m.m[2][1] * m.m[3][3] + 
              m.m[0][0]  * m.m[2][3] * m.m[3][1] + 
              m.m[2][0]  * m.m[0][1] * m.m[3][3] - 
              m.m[2][0]  * m.m[0][3] * m.m[3][1] - 
              m.m[3][0] * m.m[0][1] * m.m[2][3] + 
              m.m[3][0] * m.m[0][3] * m.m[2][1];

    inv.m[3][1] = 0;
    /*inv.m[3][1] = m.m[0][0]  * m.m[2][1] * m.m[3][2] - 
              m.m[0][0]  * m.m[2][2] * m.m[3][1] - 
              m.m[2][0]  * m.m[0][1] * m.m[3][2] + 
              m.m[2][0]  * m.m[0][2] * m.m[3][1] + 
              m.m[3][0] * m.m[0][1] * m.m[2][2] - 
              m.m[3][0] * m.m[0][2] * m.m[2][1];*/

    inv.m[0][2] = m.m[0][1]  * m.m[1][2] * m.m[3][3] - 
             m.m[0][1]  * m.m[1][3] * m.m[3][2] - 
             m.m[1][1]  * m.m[0][2] * m.m[3][3] + 
             m.m[1][1]  * m.m[0][3] * m.m[3][2] + 
             m.m[3][1] * m.m[0][2] * m.m[1][3] - 
             m.m[3][1] * m.m[0][3] * m.m[1][2];

    inv.m[1][2] = -m.m[0][0]  * m.m[1][2] * m.m[3][3] + 
              m.m[0][0]  * m.m[1][3] * m.m[3][2] + 
              m.m[1][0]  * m.m[0][2] * m.m[3][3] - 
              m.m[1][0]  * m.m[0][3] * m.m[3][2] - 
              m.m[3][0] * m.m[0][2] * m.m[1][3] + 
              m.m[3][0] * m.m[0][3] * m.m[1][2];

    inv.m[2][2] = m.m[0][0]  * m.m[1][1] * m.m[3][3] - 
              m.m[0][0]  * m.m[1][3] * m.m[3][1] - 
              m.m[1][0]  * m.m[0][1] * m.m[3][3] + 
              m.m[1][0]  * m.m[0][3] * m.m[3][1] + 
              m.m[3][0] * m.m[0][1] * m.m[1][3] - 
              m.m[3][0] * m.m[0][3] * m.m[1][1];

    inv.m[3][2] = 0;
    /*inv.m[3][2] = -m.m[0][0]  * m.m[1][1] * m.m[3][2] + 
               m.m[0][0]  * m.m[1][2] * m.m[3][1] + 
               m.m[1][0]  * m.m[0][1] * m.m[3][2] - 
               m.m[1][0]  * m.m[0][2] * m.m[3][1] - 
               m.m[3][0] * m.m[0][1] * m.m[1][2] + 
               m.m[3][0] * m.m[0][2] * m.m[1][1];*/

    inv.m[0][3] = -m.m[0][1] * m.m[1][2] * m.m[2][3] + 
              m.m[0][1] * m.m[1][3] * m.m[2][2] + 
              m.m[1][1] * m.m[0][2] * m.m[2][3] - 
              m.m[1][1] * m.m[0][3] * m.m[2][2] - 
              m.m[2][1] * m.m[0][2] * m.m[1][3] + 
              m.m[2][1] * m.m[0][3] * m.m[1][2];

    inv.m[1][3] = m.m[0][0] * m.m[1][2] * m.m[2][3] - 
             m.m[0][0] * m.m[1][3] * m.m[2][2] - 
             m.m[1][0] * m.m[0][2] * m.m[2][3] + 
             m.m[1][0] * m.m[0][3] * m.m[2][2] + 
             m.m[2][0] * m.m[0][2] * m.m[1][3] - 
             m.m[2][0] * m.m[0][3] * m.m[1][2];

    inv.m[2][3] = -m.m[0][0] * m.m[1][1] * m.m[2][3] + 
               m.m[0][0] * m.m[1][3] * m.m[2][1] + 
               m.m[1][0] * m.m[0][1] * m.m[2][3] - 
               m.m[1][0] * m.m[0][3] * m.m[2][1] - 
               m.m[2][0] * m.m[0][1] * m.m[1][3] + 
               m.m[2][0] * m.m[0][3] * m.m[1][1];

    inv.m[3][3] = 1;
    /*inv.m[3][3] = m.m[0][0] * m.m[1][1] * m.m[2][2] - 
              m.m[0][0] * m.m[1][2] * m.m[2][1] - 
              m.m[1][0] * m.m[0][1] * m.m[2][2] + 
              m.m[1][0] * m.m[0][2] * m.m[2][1] + 
              m.m[2][0] * m.m[0][1] * m.m[1][2] - 
              m.m[2][0] * m.m[0][2] * m.m[1][1];*/

    det = m.m[0][0] * inv.m[0][0] + m.m[0][1] * inv.m[1][0] 
    	+ m.m[0][2] * inv.m[2][0] + m.m[0][3] * inv.m[3][0];

    bool hasInverse;
    if (det == 0){
       hasInverse = false;
       cout << "LA MATRIZ NO ES INVERSIBLE WTF" << endl;
    }
    else{
    	det = 1.0 / det;
    	hasInverse = true;
    }
    if(hasInverse){
	    for (int i = 0; i < 4; i++){
	    	for(int j = 0; j < 4; j++){
	    		inv.m[i][j] = inv.m[i][j] * det;
	    	}
	    }
	}
    return inv;
}

} // namespace m4d

#endif // __MATRIX4X4D__
