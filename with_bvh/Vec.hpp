#ifndef Vec_hpp
#define Vec_hpp

#include <cmath>

struct Vec {
	double x, y, z;
	inline Vec(double x0=0, double y0=0, double z0=0){ x=x0; y=y0; z=z0; }
	inline Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); }
	inline Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); }
	inline Vec operator*(double b) const { return Vec(x*b,y*b,z*b); }
	inline Vec operator/(double b) const { return Vec(x/b,y/b,z/b); }
	inline Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); }
	inline Vec& norm() { return *this = *this * (1/std::sqrt(x*x+y*y+z*z)); }
	inline double length() { return std::sqrt(x*x+y*y+z*z); }
	inline double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; }
	inline Vec operator%(const Vec &b) const {return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
//	inline double& operator[](size_t i) { return data[i]; }
	inline const double& operator[](size_t i) const { return i==0 ? x : (i==1 ? y : z); }
};

#endif
