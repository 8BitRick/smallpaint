#include <string>
#include <omp.h>
#include <ctime>
#include <unordered_map>
#include <random>

using namespace std;

// Surface render settings
#define SURFACES			// render surfaces (if deactivated renders a single point light in the center of the scene and without the cornell box etc.)
#define CORNELLBOX			// render cornellbox (only if surfaces are enabled. if enabled renders a cornelbox for the scene)
// Media render settings
#define MULTIPLE_SCATTERING	// use multiple scattering
#define RM					// if enabled uses ray marching for heterogeneous media, if disabled uses delta tracking
//#define HOM					// if enabled ignores RM command and computes exponential homogeneous transmittances


// Random Number Generation
std::mt19937 mersenneTwister;
std::uniform_real_distribution<double> uniform;
#define RND (2.0*uniform(mersenneTwister)-1.0)	// [-1, 1]
#define RND2 (uniform(mersenneTwister))			// [ 0, 1]

// Constants
#define PI 3.1415926536
#define INV_4PI 0.0795774715	// 1/(4pi)
const double inf = 1e9;
const double eps = 1e-6;


// Assertions for debugging
#define ASSERT
void assert(const bool& condition, string message) {
#ifndef ASSERT
	return;
#endif
	if (!condition)
		printf((message + "\n").c_str());
}

// Basic Vector class
struct Vec {
	double x, y, z;
	Vec(double x0, double y0, double z0){ x = x0; y = y0; z = z0; }
	Vec(double xyz0 = 0){ x = xyz0; y = xyz0; z = xyz0; }
	Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
	Vec operator+=(const Vec &b) { x += b.x; y += b.y; z += b.z; return (*this); }
	Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
	Vec operator*(double b) const { return Vec(x*b, y*b, z*b); }
	Vec operator*(const Vec &b) const { return Vec(x*b.x, y*b.y, z*b.z); }
	Vec operator*=(const Vec &b) { x *= b.x; y *= b.y; z *= b.z; return (*this); }
	Vec operator/(double b) const { return Vec(x / b, y / b, z / b); }
	Vec operator/(const Vec &b) const { return Vec(x / b.x, y / b.y, z / b.z); }
	bool operator<(const Vec &b) const { return x < b.x && y < b.y && z < b.z; }
	bool operator>(const Vec &b) const { return x > b.x && y > b.y && z > b.z; }
	Vec& norm(){ return *this = *this * (1 / sqrt(x*x + y*y + z*z)); }
	double length() const { return sqrt(x*x + y*y + z*z); }
	double dot(const Vec &b) const { return x*b.x + y*b.y + z*b.z; }
	double avg() const { return (x + y + z) / 3.0; }
	double max() const { return x > y ? (x > z ? x : z) : (y > z ? y : z); }
	double min() const { return x < y ? (x < z ? x : z) : (y < z ? y : z); }
	Vec operator%(const Vec &b) const { return Vec(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x); }
	const double& operator[](size_t i) const { return i == 0 ? x : (i == 1 ? y : z); }
};
Vec operator*(double a, const Vec &b) { return Vec(a * b.x, a * b.y, a * b.z); } // double * Vec


// Helper functions
#pragma region helpers

// Create an orthonormal system
// given v1, set v2 and v3 so they form an orthonormal system
// (we assume v1 is already normalized)
void ons(const Vec& v1, Vec& v2, Vec& v3) {
	if (std::abs(v1.x) > std::abs(v1.y)) {
		// project to the y = 0 plane and construct a normalized orthogonal vector in this plane
		float invLen = 1.f / sqrtf(v1.x * v1.x + v1.z * v1.z);
		v2 = Vec(-v1.z * invLen, 0.0f, v1.x * invLen);
	}
	else {
		// project to the x = 0 plane and construct a normalized orthogonal vector in this plane
		float invLen = 1.0f / sqrtf(v1.y * v1.y + v1.z * v1.z);
		v2 = Vec(0.0f, v1.z * invLen, -v1.y * invLen);
	}
	v3 = v1 % v2;
}

// Linear interpolation
inline double lerp(const double& a, const double& b, const double alpha) {
	return (1.0 - alpha) * a + alpha * b;
}

// Convert double [0,1] to int [0,255] with possible gamma correction
int clrToInt(double clr, bool gammacorrection = false) {
	if (gammacorrection) clr = pow(clr, 1.0 / 2.2);
	int converted = (int)(clr * 255 + 0.5);
	return converted < 0 ? 0 : converted > 255 ? 255 : converted;
}

#pragma endregion

// Perlin noise
#pragma region perlin

// Random noise value by 3D-seed
double noise(int x, int y, int z)
{
	int n = x + y * 57 + z * 101;
	n = (n << 13) ^ n;
	int t = (n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff;
	return 1.0 - double(t) * 0.731322574615478515625e-9;
}

// Cosine interpolation
double cosintp(double a, double b, double alpha)
{
	double x = (1.0 - cos(alpha * PI)) * 0.5;
	return a * (1 - x) + b * x;
}

// Trilinear noise interpolation
double intpnoise(double x, double y, double z)
{

	unsigned int x0 = (int)floor(x);
	unsigned int y0 = (int)floor(y);
	unsigned int z0 = (int)floor(z);
	unsigned int x1 = (int)ceil(x);
	unsigned int y1 = (int)ceil(y);
	unsigned int z1 = (int)ceil(z);
	double xd = (x - x0);
	double yd = (y - y0);
	double zd = (z - z0);

	// cosine interpolation along x
	double c00 = cosintp(noise(x0, y0, z0), noise(x1, y0, z0), xd);
	double c10 = cosintp(noise(x0, y1, z0), noise(x1, y1, z0), xd);
	double c01 = cosintp(noise(x0, y0, z1), noise(x1, y0, z1), xd);
	double c11 = cosintp(noise(x0, y1, z1), noise(x1, y1, z1), xd);

	// along y
	double c0 = cosintp(c00, c10, yd);
	double c1 = cosintp(c01, c11, yd);

	// along z
	double c = cosintp(c0, c1, zd);

	return c;
}

// Perlin noise at specified position
// Output range is approximately -1 to 1
double perlin(double x, double y, double z)
{
	int octaves = 8;
	double persistence = 0.25;
	double total = 0.0;
	for (int i = 0; i < octaves; ++i)
	{
		double freq = pow(2, i);
		double ampl = pow(persistence, i);
		total += intpnoise(x * freq, y * freq, z * freq) * ampl;
	}

	return total;
}

#pragma endregion

// Grid
#pragma region grid
// 3D-Grid of double values
class Grid {

public:
	unsigned int n; // Number of elements in one dimension --> Grid is n by n by n
	double* values;

	Grid(unsigned int _n = 0, double* _v = 0)
	{
		n = _n;
		values = _v;
	}

	double& el(unsigned int x, unsigned int y, unsigned int z)
	{
		int index = z * n * n + y * n + x;
		if (index > n*n*n - 1)
			printf("Grid index out of range");
		return values[index];
	}

	double& operator()(unsigned x, unsigned y, unsigned z) { return el(x, y, z); }
};

// Method for filling a grid with values (uniform random, 2D checkerboard, x gradient and 3D perlin noise)
void fillGrid(const int type, Grid& grid)
{
	const unsigned int n = grid.n;

	if (type == 0) // uniform random
	{
		for (int i = 0; i < n*n*n; ++i)
			grid.values[i] = RND2;
	}
	else if (type == 1) // xy checkerboard
	{
		for (int x = 0; x < n; ++x)
		{
			for (int y = 0; y < n; ++y)
			{
				for (int z = 0; z < n; ++z)
				{
					grid(x, y, z) = ((x + y) % 2 == 1) ? 1.0 : 0.0;
				}
			}
		}
	}
	else if (type == 2) // x gradient
	{
		for (int x = 0; x < n; ++x)
		{
			for (int y = 0; y < n; ++y)
			{
				for (int z = 0; z < n; ++z)
				{
					grid(x, y, z) = sqrt(x*x) / sqrt(n*n);
				}
			}
		}
	}
	else if (type == 3) // 3d perlin noise
	{
		for (int x = 0; x < n; ++x)
		{
			for (int y = 0; y < n; ++y)
			{
				for (int z = 0; z < n; ++z)
				{
					double p = perlin(x, y, z); //more or less -1,1
					double pp = (p + 1.0) * 0.5; // map to 0,1
					double ppp = pp > 1.0 ? 1.0 : pp < 0.0 ? 0.0 : pp; // clamp to 0,1
					grid(x, y, z) = ppp;
				}
			}
		}
	}
	else if (type == 4) // homogeneous
	{
		for (int i = 0; i < n*n*n; ++i)
			grid.values[i] = 1;
	}
}
#pragma endregion

// Scene structs/classes/methods
#pragma region scene
// Ray with origin and (normalized) direction
struct Ray {
	Vec o, d;
	Ray(Vec o0 = 0, Vec d0 = 0) { o = o0, d = d0.norm(); }
};

// Sceneobject with a material
class Obj {
public:
	Vec clr;
	double emission;
	int type; // 1 = diffuse, 2 = specular, 3 = refractive

	void setMat(Vec clr_ = 0, double emission_ = 0, int type_ = 0) { clr = clr_; emission = emission_; type = type_; }

	// Subclasses have to implement an intersection routine and normal computation
	virtual double intersect(const Ray&) const = 0;
	virtual Vec normal(const Vec&) const = 0;
};

// Plane defined by normal and distance
class Plane : public Obj {
public:
	Vec n;
	double d;
	Plane(double d_ = 0, Vec n_ = 0) {
		d = d_;
		n = n_;
	}
	double intersect(const Ray& ray) const {
		double d0 = n.dot(ray.d);
		if (d0 != 0) {
			double t = -1 * (((n.dot(ray.o)) + d) / d0);
			return (t > eps) ? t : 0;
		}
		else return 0;
	}
	Vec normal(const Vec& p0) const { return n; }
};

// Sphere defined by center and radius
class Sphere : public Obj {
public:
	Vec c;
	double r;

	Sphere(double r_ = 0, Vec c_ = 0) { c = c_; r = r_; }
	double intersect(const Ray& ray) const { //Note: assuming normalized direction
		double b = ((ray.o - c) * 2).dot(ray.d);
		double c_ = (ray.o - c).dot((ray.o - c)) - (r*r);
		double disc = b*b - 4 * c_;
		if (disc<0) return 0;
		else disc = sqrt(disc);
		double sol1 = -b + disc;
		double sol2 = -b - disc;
		return (sol2>eps) ? sol2 / 2 : ((sol1>eps) ? sol1 / 2 : 0);
	}

	Vec normal(const Vec& p0) const {
		return (p0 - c).norm();
	}
};

// Intersection class storing sceneobject reference and ray distance parameter t
class Intersection {
public:
	Intersection() { t = inf; object = nullptr; }
	Intersection(double t_, Obj* object_) { t = t_; object = object_; }
	operator bool() { return object != nullptr; }
	double t;
	Obj* object;
};

// Pointlight with a color, intensity and position
class PointLight {

private:
	Vec clr;
	Vec pos;
	double intensity;

public:
	PointLight(Vec pos_ = 0, double intensity_ = 0, Vec clr_ = 0) {
		pos = pos_;
		intensity = intensity_;
		clr = clr_;
	}

	Vec emission() {
		return clr * intensity;
	}

	Vec samplePos()	{
		return pos;
	}
};

// Heterogeneous medium
// The medium is an n by n by n grid of densities with positions parametrized to [0,1]^3
// I.e. densities[0][0][0] = left bottom front, densities[1][1][1] = right top back
class Medium {

private:
	Vec _sigma_a;	// base absorption coefficient
	Vec _sigma_s;	// base scattering coefficient
	Vec _sigma_t;	// base extinction coefficient
	Vec _emission;	// base volumetric emission

public:
	Vec position;	// left bottom front corner position
	Vec size;		// size of the medium
	Grid densities;	// density value grid

	Medium(Vec _a = 0, Vec _s = 0, Vec _em = 0, Grid _densities = 0, Vec _pos = 0, Vec _size = 0)	{
		_sigma_a = _a;
		_sigma_s = _s;
		_sigma_t = _a + _s;
		_emission = _em;
		densities = _densities;
		position = _pos;
		size = _size;
	}

	// get the density at a position p
	double density(Vec p)
	{
		if (densities.n < 2) return 1.0; //grid too small

		//p is in worldspace, medium in worldspace is in [position,position+size]
		// -> p in param space: (p-position) / size ... [0, 1]
		// -> to index space: * n ... [0, n-1]
		Vec p_is = (p - position) / size * (densities.n - 1); // p in index space [0, n-1]^3

		//Note: assuming that position p is inside the medium --> no out of bounds check

		//Find nearest grid values and interpolate

		unsigned int x0 = (int)floor(p_is.x); //prev neighbour index
		unsigned int y0 = (int)floor(p_is.y);
		unsigned int z0 = (int)floor(p_is.z);
		unsigned int x1 = (int)ceil(p_is.x); //next neighbour index
		unsigned int y1 = (int)ceil(p_is.y);
		unsigned int z1 = (int)ceil(p_is.z);
		double xd = (p_is.x - x0); //Note: x1 - x0 == 1, division by (x1 - x0) omitted
		double yd = (p_is.y - y0);
		double zd = (p_is.z - z0);

		// lerp along x
		double c00 = lerp(densities(x0, y0, z0), densities(x1, y0, z0), xd);
		double c10 = lerp(densities(x0, y1, z0), densities(x1, y1, z0), xd);
		double c01 = lerp(densities(x0, y0, z1), densities(x1, y0, z1), xd);
		double c11 = lerp(densities(x0, y1, z1), densities(x1, y1, z1), xd);

		// lerp along y
		double c0 = lerp(c00, c10, yd);
		double c1 = lerp(c01, c11, yd);

		//lerp along z
		double c = lerp(c0, c1, zd);

		return c;
	}

	// Medium coefficients at a position p
	Vec sigma_a(Vec p) { return _sigma_a * density(p); }
	Vec sigma_s(Vec p) { return _sigma_s * density(p); }
	Vec sigma_t(Vec p) { return _sigma_t * density(p); }
	Vec emission(Vec p) { return _emission * density(p); }

	// Check if point p lies within the medium
	bool inbounds(const Vec& p)
	{
		Vec a = position;
		Vec b = position + size;
		return (p.x >= a.x && p.x <= b.x && p.y >= a.y && p.y <= b.y && p.z >= a.z && p.z <= b.z);
	}

	// Clamp starting and ending point t0 and t1 of a ray to the medium
	void clampRay(const Ray& ray, double& t0, double& t1)
	{
		//Model medium sides with 6 planes
		//intersect ray with planes
		//map t0 to nearest intersection (if there are two, else t0 lies inside the medium already)
		//map t1 to farthest intersection
		//set t0 to t1 if no intersection

		Vec b = position + size; // medium end

		double tsides[6];
		tsides[0] = Plane(-position.x, Vec(1, 0, 0)).intersect(ray);
		tsides[1] = Plane(b.x, Vec(-1, 0, 0)).intersect(ray);
		tsides[2] = Plane(-position.y, Vec(0, 1, 0)).intersect(ray);
		tsides[3] = Plane(b.y, Vec(0, -1, 0)).intersect(ray);
		tsides[4] = Plane(-position.z, Vec(0, 0, 1)).intersect(ray);
		tsides[5] = Plane(b.z, Vec(0, 0, -1)).intersect(ray);

		//limit planes to grid/cube-faces
		for (int i = 0; i < 6; ++i)
		{
			//point on the plane
			Vec pi = ray.o + tsides[i] * ray.d;
			//any plane intersection producing a point outside of the grid-cube gets t set to 0
			if (pi.x < position.x || pi.y < position.y || pi.z < position.z
				|| pi.x > b.x || pi.y > b.y || pi.z > b.z)
				tsides[i] = 0;
		}

		// find nearest and farthest intersection
		double t_nearest = inf, t_farthest = 0;
		for (int i = 0; i < 6; ++i)
		{
			if (tsides[i] > 0 && tsides[i] < inf)
			{
				if (tsides[i] < t_nearest) t_nearest = tsides[i];
				if (tsides[i] > t_farthest) t_farthest = tsides[i];
			}
		}

		// remap t0 and t1 if necessary
		if (t_nearest < inf && t_farthest > 0) //intersection
		{
			if (t_nearest == t_farthest) //only one intersection (only exit)
			{
				if (t_farthest < t1)
					t1 = t_farthest;
			}
			else // two intersections (entry and exit)
			{
				if (t_farthest < t1)
					t1 = t_farthest;
				if (t_nearest > t0)
					t0 = t_nearest;
			}
		}
		else
		{
			t0 = t1; //no intersection
		}

	}

	// Transmittance inside a homogeneous unbounded medium for a given travelled distance
	Vec transmittanceHom(double distance) {
		Vec tau = _sigma_t * distance;
		return Vec(exp(-tau.x), exp(-tau.y), exp(-tau.z));
	}

	// delta/woodcock tracking transmittance
	Vec transmittanceDelta(const Ray& ray, const double &tmin, const double &tmax, const int num_samples = 2) {

		//Note: delta tracking computes transmittance as a diraq function, i.e. either completely 1 or 0.
		// --> by iteratively applying the method for a number of samples and averaging the result the real transmittance is estimated

		if (num_samples <= 0) return Vec(1.0);

		//clamp the ray to the bounded medium 
		double tmediamin = tmin, tmediamax = tmax;
		clampRay(ray, tmediamin, tmediamax);

		if (tmediamin >= tmediamax) //medium not intersected
			return Vec(1.0);

		//averaged result
		double result = 0;

		//maximum/majorant density of the medium along the ray
		double max_density = majorantDensity(ray, tmediamin, tmediamax);

		if (max_density == 0) //ray with max density = 0 --> full transmittance
			return Vec(1.0);

		const double inv_max_density = 1 / max_density;
		//Note: sigma_t is an rgb vector but we need a 1D value.
		// We found no discussion in the respective papers as to how to work with rgb data
		// Therefore we just use the maximum color. Other options would be, e.g., to use an rgb-to-luminance conversion
		const double inv_max_sigma_t = inv_max_density / _sigma_t.max();

		for (int i = 0; i < num_samples; ++i) {
			double t = tmediamin; //reset  the position  to the beginning of the ray or medium
			while (true) {
				t -= log(1 - RND2) * inv_max_sigma_t; //randomly step through the medium assuming it has the majorant density everywhere
				if (t >= tmediamax) { //the end was reached and no real particle collision happened
					result += 1.0; //--> full transmittance
					break;
				}

				double d = density(ray.o + t * ray.d);
				if (RND2 <= d * inv_max_density) //check if a real particle collision happened. Note: sigma_t can be omitted here because it is canceled out in multiplication by the inverse anyway
					break; //-->zero transmittance, because a "real" particle was hit
			}
		}

		// return result averaged over number of samples
		return Vec(result / num_samples);
	}

	// Transmittance in a heterogeneous medium by ray marching
	Vec transmittanceRM(const Ray& ray, double t0, double t1) {

		// Initial tau
		Vec tau(0.0);

		// find starting and ending point of the ray in the medium
		clampRay(ray, t0, t1);

		if (t0 >= t1) // no medium was intersected
			return Vec(1);

		// Calculate a stepsize
		double stepsize = 0.1;
		double numsteps = ceil((t1 - t0) / stepsize);
		if (numsteps == 0) return 1.0; // no steps taken --> full transmittance
		double step = (t1 - t0) / numsteps;

		for (double t = t0/* + RND2 * 0.9 * step*/; t <= t1; t += step) // Note: possible random offseting
		{
			// Calculate tau for the current stepping position and add to accumulated tau
			Vec tau_step = sigma_t(ray.o + t * ray.d);
			tau += tau_step;

			// Note: possible russian roulette termination when transmittance low
		}

		tau *= step; //weigh by step size

		// calculate exponential transmittance using the accumulated tau and return it
		return Vec(exp(-tau.x), exp(-tau.y), exp(-tau.z));
	}

	// calculate the majorant/maximum density along a ray
	double majorantDensity(const Ray& ray, const double& t0, const double& t1)
	{
		//majorant density
		double m = 0;

		//Note: Assuming t0 and t1 are clamped to medium bounds

		if (t0 >= t1)
			return 0; //Medium not intersected

		double stepsize = 0.1;
		double length = (t1 - t0);
		double numsteps = ceil(length / stepsize);
		if (numsteps == 0)
			numsteps = 1;
		double step = (t1 - t0) / numsteps;

		for (double t = t0/* + RND2 * 0.9 * step*/; t <= t1; t += step) //Note: possible random offseting
		{
			// update maximum density
			double tmp = density(ray.o + t * ray.d);
			if (tmp > m)
				m = tmp;
		}

		return m;
	}
};

class Scene {
	vector<Obj*> objects;

public:
	PointLight* light;
	Medium* medium;

	void add(Obj* object) {
		objects.push_back(object);
	}

	//naively intersect all objects in the scene and return the closest intersection
	Intersection intersect(const Ray& ray) const {
		Intersection closestIntersection;
		// intersect all objects, one after the other
		for (auto iter = objects.begin(); iter != objects.end(); ++iter) {
			double t = (*iter)->intersect(ray);
			if (t > eps && t < closestIntersection.t) {
				closestIntersection.t = t;
				closestIntersection.object = *iter;
			}
		}
		return closestIntersection;
	}

};

#pragma endregion

// Sampling methods
#pragma region sampling

// Uniform sampling on a hemisphere to produce outgoing ray directions.
// courtesy of http://www.rorydriscoll.com/2009/01/07/better-sampling/
Vec sampleHemisphere(double u1, double u2) {
	const double r = sqrt(1.0 - u1*u1);
	const double phi = 2 * PI * u2;
	return Vec(cos(phi)*r, sin(phi)*r, u1);
}

// Uniform sampling on a complete unit sphere
Vec sampleSphere(double u1, double u2) {
	const double u1d = u1 * 2.0 - 1.0; // remap to -1 .. 1
	const double r = sqrt(1.0 - u1d*u1d);
	const double phi = 2 * PI * u2;
	return Vec(cos(phi)*r, sin(phi)*r, u1d);
}

// Uniform sampling along a ray
void sampleUniform(double u, double maxDistance, double &distance, double& probability) {
	distance = u * maxDistance;
	probability = 1.0 / maxDistance;
}

// Exponential sampling along ray
void sampleExponential(double u, double maxDistance, const Vec& sigma_t, double &distance, double& probability) {

	double sigma = (sigma_t.x + sigma_t.y + sigma_t.z) / 3.0; // using average

	// remap u to account for finite max distance
	double minU = exp(-sigma*maxDistance);
	double a = u*(1.0 - minU) + minU;

	// sample with pdf proportional to exp(-sig*d) // see http://www.cs.cornell.edu/courses/cs6630/2012sp/notes/09volpath.pdf
	distance = -log(a) / sigma;
	probability = sigma*a / (1.0 - minU);
}

// Equi-angular sampling along a ray w.r.t. the light source
void sampleEquiAngular(double u, double maxDistance, const Ray& ray, const Vec& lightPos, double &distance, double& probability) {
	// Get closest point to light along the ray
	double tclose = (lightPos - ray.o).dot(ray.d);

	// Closest distance to light from ray
	double rayDist = (ray.o + ray.d * tclose - lightPos).length();

	// angles spanned by light, closest point and ray origin or end
	double thetaA = atan((0.0 - tclose) / rayDist);
	double thetaB = atan((maxDistance - tclose) / rayDist);

	double phi = (1 - u) * thetaA + u * thetaB;

	// sample ray position based on angle
	double t = rayDist * tan(phi); // relative to closest point
	distance = tclose + t;
	probability = rayDist / ((thetaB - thetaA)*(rayDist*rayDist + t*t)); //see: https://www.solidangle.com/research/s2011_equiangular_slides.pdf
}

// Sample image plane by pixel coordinates
// Input is the pixel offset, output is the appropriate coordinate
// on the image plane
Vec camcr(const double x, const double y, const double width, const double height) {
	double w = width;
	double h = height;
	float fovx = PI / 4;
	float fovy = (h / w) * fovx;
	return Vec(((2 * x - w) / w) * tan(fovx),
		-((2 * y - h) / h) * tan(fovy),
		-1.0);
}

#pragma endregion
