#include <Matrix.h>
#include <Point.h>
#include <Vector.h>
#include <Geometry_Common.h>
#include <iomanip>
#include <cmath>

namespace geometry{

double Matrix::Det3::operator()() const{
	return data[0][0] * data[1][1] * data[2][2] + data[0][1] * data[1][2] * data[2][0] + data[0][2] * data[1][0] * data[2][1] -
		data[0][0] * data[1][2] * data[2][1] - data[0][1] * data[1][0] * data[2][2] - data[0][2] * data[1][1] * data[2][0] ;
}

Matrix::Matrix(){
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++){
			data[i][j] = 0 ;
		}
}

Matrix::Matrix(const double d00, const double d01, const double d02, const double d03,
	       const double d10, const double d11, const double d12, const double d13,
	       const double d20, const double d21, const double d22, const double d23,
	       const double d30, const double d31, const double d32, const double d33){
	set(d00, d01, d02, d03, d10, d11, d12, d13, d20, d21, d22, d23, d30, d31, d32, d33) ;
}

void Matrix::set(const double d00, const double d01, const double d02, const double d03,
		 const double d10, const double d11, const double d12, const double d13,
		 const double d20, const double d21, const double d22, const double d23,
		 const double d30, const double d31, const double d32, const double d33){
	data[0][0] = d00, data[0][1] = d01, data[0][2] = d02, data[0][3] = d03 ;
	data[1][0] = d10, data[1][1] = d11, data[1][2] = d12, data[1][3] = d13 ;
	data[2][0] = d20, data[2][1] = d21, data[2][2] = d22, data[2][3] = d23 ;
	data[3][0] = d30, data[3][1] = d31, data[3][2] = d32, data[3][3] = d33 ;
}

double *Matrix::operator[](const int rowNum){
	if (rowNum < 0 || rowNum >= 4){
		throw IndexOutOfRange("Matrix::operator[]") ;
	}
	return data[rowNum] ;
}

const double *Matrix::operator[](const int rowNum) const{
	if (rowNum < 0 || rowNum >= 4){
		throw IndexOutOfRange("Matrix::operator[]") ;
	}
	return data[rowNum] ;
}

bool Matrix::isZero() const{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++){
			if (data[i][j] != 0)
				return false ;
		}
	return true ;
}

bool Matrix::isIdentity() const{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++){
			if (i == j){
				if (data[i][j] != 1)
					return false ;
			}
			else{
				if (data[i][j] != 0)
					return false ;
			}
		}
	return true ;
}

double Matrix::determinant() const{
	double temp = 0 ;
	for (int i = 0; i < 4; i++){
		Det3 det3 = getDet3(0, i) ;
		temp += ((i % 2 == 0)? 1 : -1) * data[0][i] * det3() ;
	}
	return temp ;
}

Matrix::Det3 Matrix::getDet3(const int r, const int c) const{
	if (r < 0 || r >= 4 || c < 0 || c >= 4){
		throw IndexOutOfRange("Matrix::getDet3") ;
	}
	Det3 det3 ;
	int rn = -1 ;
	for (int i = 0; i < 4 ;i++){
		if (i == r)
			continue ;
		rn++ ;
		int cn = -1 ;
		for (int j = 0; j < 4; j++){
			if (j == c)
				continue ;
			cn++ ;
			det3.data[rn][cn] = data[i][j] ;
		}
	}
	return det3 ;
}

std::ostream &operator<<(std::ostream &os, const Matrix &m){
	for (int i = 0; i < 4; i++){
		if (i == 0){
			os << "Matrix: [" << std::fixed ;
		}
		else{
			os << std::setw(9) << ' ' ;
		}
		for (int j = 0; j < 4; j++){
			os << std::setw(13) << m[i][j] ;
		}
		if (i == 3){
			os << ']' ;
		}
		else{
			os << std::endl ;
		}
	}
	return os ;
}

bool operator==(const Matrix &m1, const Matrix &m2){
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			if (m1[i][j] != m2[i][j])
				return false ;
	return true ;
}

Matrix operator*(const Matrix &m1, const Matrix &m2){
	Matrix r ;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++){
			double temp = 0 ;
			for (int k = 0; k < 4; k++){
				temp += m1[i][k] * m2[k][j] ;
			}
			r[i][j] = temp ;
		}
	return r ;
}

Matrix inverse(const Matrix& m){
	Matrix r ;
	double matrixDeterminant = m.determinant() ;
	if (matrixDeterminant == 0)	//The determinant is zero, means that there does not exist an inverse matrix
		return r ;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++){
			Matrix::Det3 det3 = m.getDet3(i, j) ;
			r[j][i] = (((i + j) % 2 == 0)? 1 : -1) * det3() / matrixDeterminant ;
		}
	return r ;
}

Matrix identity(){
	Matrix r ;
	for (int i = 0; i < 4; i++)
		r[i][i] = 1 ;
	return r ;
}

Matrix lookAt(const Point &cameraPos ,const Point &centerPos, const Vector &up){
	Vector view(cameraPos, centerPos) ;
	view.normalize() ;		//Crucial: normalize the vectors to make them a frame
	Vector backVector = -view ;
	Vector upVector(up) ;
	upVector.normalize() ;
	Vector rightVector = crossProduct(view, upVector) ;
	Matrix m ;
	for (int i = 0; i < 3; i++){
		m[i][0] = rightVector[i] ;
		m[i][1] = upVector[i] ;
		m[i][2] = backVector[i] ;
		m[i][3] = cameraPos[i] ;
	}
	m[3][3] = 1 ;
	return inverse(m) ;
}

Matrix orthoProjection(const double left, const double right, const double bottom, const double top, const double near, const double far){
	//NOTICE that near and far are DISTANCES to the near and far depth clipping planes, so if they are both positive, near < far!!
	//Since we want to put camera at origin, pointing to -z, the actual clipping planes are -near and -far!
	//The Ortho Projection matrix is S * T (left-multiply matrix), and S means scale, T means translate
	//T matrix: move center of the view volume to the ORIGIN, {-(left + right)/2, -(bottom + top)/2, (near + far)/2}
	//NOTE that since the actual clip planes are -near and -far, we are moving the view volume along +z! 
	//S matrix: scale to make the view volume have sides of length 2, {2/(right - left), 2/(top - bottom), -2/(far - near)}
	//NOTE that an additional minus sign is used when scaling z, since after the translation, the near clipping plane
	//will be on the positive side of the z axis, while the far clipping plane will be on the negative side of the z
	//axis. This is inconvenient for depth computation because nearer points have z components larger than farther points
	//do. So, we add a minus sign to reverse it, make nearer points have smaller z values. NOTE that projection is the
	//last step before outputting pixels, so this inversion step is a must!
	//This is to normalize the projection, namely put the view volume in cube {(-1, -1, -1), (1, 1, 1)}
	return Matrix(
		2 / (right - left),	0,			0,			-(right + left) / (right - left),
		0,			2 / (top - bottom),	0,			-(top + bottom) / (top - bottom),
		0,			0,			-2 / (far - near),	-(far + near) / (far - near),
		0,			0,			0,			1
	) ;
}

Matrix rotateAroundX(const double angle){
	double theta = angle / 180 * PI ;
	return Matrix(
		1,			0,			0,			0,
		0,			std::cos(theta),	-std::sin(theta),	0,
		0,			std::sin(theta),	std::cos(theta),	0,
		0,			0,			0,			1
	) ;
}

Matrix rotateAroundY(const double angle){
	double theta = angle / 180 * PI ;
	return Matrix(
		std::cos(theta),	0,			std::sin(theta),	0,
		0,			1,			0,			0,
		-std::sin(theta),	0,			std::cos(theta),	0,
		0,			0,			0,			1
	) ;
}

Matrix rotateAroundZ(const double angle){
	double theta = angle / 180 * PI ;
	return Matrix(
		std::cos(theta),	-std::sin(theta),	0,			0,
		std::sin(theta),	std::cos(theta),	0,			0,
		0,			0,			1,			0,
		0,			0,			0,			1
	) ;
}

Matrix translate(const double tx, const double ty, const double tz){
	return Matrix(
		1,		0,		0,		tx,
		0,		1,		0,		ty,
		0,		0,		1,		tz,
		0,		0,		0,		1
	) ;
}

Matrix scale(const double fx, const double fy, const double fz){
	return Matrix(
		fx,		0,		0,		0,
		0,		fy,		0,		0,
		0,		0,		fz,		0,
		0,		0,		0,		1
	) ;
}

}	//End namespace geometry
