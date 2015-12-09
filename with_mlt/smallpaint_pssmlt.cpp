// smallpaint by karoly zsolnai - zsolnai@cg.tuwien.ac.at
//
// compilation by: g++ smallpaint.cpp -O3 -std=gnu++0x -fopenmp
// uses >=gcc-4.5.0
// render, modify, create new scenes, tinker around, and most of all:
// have fun!
//
// This program is used as an educational learning tool on the Rendering
// course at TU Wien. Course webpage:
// http://cg.tuwien.ac.at/courses/Rendering/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <ctime>
#include <vector>
#include <string>
#include <unordered_map>
#include <random>
#include <cstdint>
#include <algorithm>
#include <stack>
#include <vector>

// Helpers for random number generation
std::mt19937 mersenneTwister;
std::uniform_real_distribution<double> uniform;

#define RND (2.0*uniform(mersenneTwister)-1.0)
#define RND2 (uniform(mersenneTwister))

#define PI 3.1415926536
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const int width=900, height=900;
const double inf=1e9;
const double eps=1e-6;
const float plarge=0.5f; // Probability of large step for PSSMLT
using namespace std;
typedef unordered_map<string, double> pl;

struct Vec {
	double x, y, z;
	Vec(double x0=0, double y0=0, double z0=0){ x=x0; y=y0; z=z0; }
	Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); }
	Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); }
	Vec operator*(double b) const { return Vec(x*b,y*b,z*b); }
	Vec operator/(double b) const { return Vec(x/b,y/b,z/b); }
	Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); }
	Vec& norm(){ return *this = *this * (1/sqrt(x*x+y*y+z*z)); }
	double length() { return sqrt(x*x+y*y+z*z); }
	double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; }
	Vec operator%(const Vec &b) const {return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
//	double& operator[](size_t i) { return data[i]; }
	const double& operator[](size_t i) const { return i==0 ? x : (i==1 ? y : z); }
};

// create an orthonormal system, assuming v1 is already normalized
void ons(const Vec& v1, Vec& v2, Vec& v3) {
    if (std::abs(v1.x) > std::abs(v1.y)) {
		// project to the y = 0 plane and construct a normalized orthogonal vector in this plane
		float invLen = 1.f / sqrtf(v1.x * v1.x + v1.z * v1.z);
		v2 = Vec(-v1.z * invLen, 0.0f, v1.x * invLen);
    } else {
		// project to the x = 0 plane and construct a normalized orthogonal vector in this plane
		float invLen = 1.0f / sqrtf(v1.y * v1.y + v1.z * v1.z);
		v2 = Vec(0.0f, v1.z * invLen, -v1.y * invLen);
    }
    v3 = v1 % v2;
}

// Rays have origin and direction.
// The direction vector should always be normalized.
struct Ray {
	Vec o, d;
	Ray(Vec o0 = 0, Vec d0 = 0) { o = o0, d = d0.norm(); }
};

struct AABB {
	AABB() { min = Vec(inf, inf, inf); max = Vec(-inf, -inf, -inf); }	// an empty interval
	AABB(Vec min_, Vec max_) { min = min_; max = max_; }
	bool unbounded() const { return min.x == -inf || min.y == -inf || min.z == -inf || max.x == inf || max.y == inf || max.z == inf; }
	size_t largestDimension() const {
		double dx = abs(max.x - min.x);
		double dy = abs(max.y - min.y);
		double dz = abs(max.z - min.z);
		if (dx > dy && dx > dz) {
			return 0;
		}
		if (dy > dz) {
			return 1;
		}
		return 2;
	}

	// ray-slab tests, see PBRT 2nd edition, section 4.2.1
	bool intersect(const Ray& ray, const Vec& inverseDirection, double closestKnownT) const {
		bool xDirNegative = ray.d.x < 0;
		bool yDirNegative = ray.d.y < 0;
		bool zDirNegative = ray.d.z < 0;

		// check for ray intersection against x and y slabs
		float tmin = ((xDirNegative ? max.x : min.x) - ray.o.x) * inverseDirection.x;
		float tmax = ((xDirNegative ? min.x : max.x) - ray.o.x) * inverseDirection.x;
		float tymin = ((yDirNegative ? max.y : min.y) - ray.o.y) * inverseDirection.y;
		float tymax = ((yDirNegative ? min.y : max.y) - ray.o.y) * inverseDirection.y;
		if (tmin > tymax || tymin > tmax) {
		    return false;
		}
		if (tymin > tmin) {
			tmin = tymin;
		}
		if (tymax < tmax) {
			tmax = tymax;
		}

		// check for ray intersection against z slab
		float tzmin = ((zDirNegative ? max.z : min.z) - ray.o.z) * inverseDirection.z;
		float tzmax = ((zDirNegative ? min.z : max.z) - ray.o.z) * inverseDirection.z;
		if (tmin > tzmax || tzmin > tmax) {
		    return false;
		}
		if (tzmin > tmin) {
		    tmin = tzmin;
		}
		if (tzmax < tmax) {
		    tmax = tzmax;
		}
		return (tmin < closestKnownT) && (tmax > eps);
	}

	Vec min;
	Vec max;
};

AABB enclose(const AABB& firstBoundingBox, const AABB& secondBoundingBox) {
	AABB ret;

	ret.min.x = std::min(firstBoundingBox.min.x, secondBoundingBox.min.x);
	ret.min.y = std::min(firstBoundingBox.min.y, secondBoundingBox.min.y);
	ret.min.z = std::min(firstBoundingBox.min.z, secondBoundingBox.min.z);

	ret.max.x = std::max(firstBoundingBox.max.x, secondBoundingBox.max.x);
	ret.max.y = std::max(firstBoundingBox.max.y, secondBoundingBox.max.y);
	ret.max.z = std::max(firstBoundingBox.max.z, secondBoundingBox.max.z);

	return ret;
}

AABB enclose(const AABB& boundingBox, const Vec& point) {
	AABB ret;

	ret.min.x = std::min(boundingBox.min.x, point.x);
	ret.min.y = std::min(boundingBox.min.y, point.y);
	ret.min.z = std::min(boundingBox.min.z, point.z);

	ret.max.x = std::max(boundingBox.max.x, point.x);
	ret.max.y = std::max(boundingBox.max.y, point.y);
	ret.max.z = std::max(boundingBox.max.z, point.z);

	return ret;
}

// Objects have color, emission, type (diffuse, specular, refractive)
// All object should be intersectable and should be able to compute their surface normals.
class Obj {
	public:
	Vec cl;
	double emission;
	int type;
	void setMat(Vec cl_ = 0, double emission_ = 0, int type_ = 0) { cl=cl_; emission=emission_; type=type_; }
	virtual double intersect(const Ray&) const = 0;
	virtual Vec normal(const Vec&) const = 0;
	virtual AABB getAABB() const = 0;
};

class Plane : public Obj {
	public:
	Vec n;
	double d;
	Plane(double d_ = 0, Vec n_= 0) {
		d=d_;
		n=n_;
	}
	double intersect(const Ray& ray) const {
		double d0 = n.dot(ray.d);
		if(d0 != 0) {
			double t = -1 * (((n.dot(ray.o))+d) / d0);
			return (t > eps) ? t : 0;
		}
		else return 0;
	}
	Vec normal(const Vec& p0) const { return n; }
	AABB getAABB() const {
		if (n.x == 0 && n.y == 0 ) return AABB(Vec(-inf, -inf, d * n.z), Vec(inf, inf, d * n.z));
		if (n.x == 0 && n.z == 0 ) return AABB(Vec(-inf, d * n.y, -inf), Vec(inf, d * n.y, inf));
		if (n.y == 0 && n.z == 0 ) return AABB(Vec(d * n.x, -inf, -inf), Vec(d * n.x, inf, inf));
		return AABB(Vec(-inf, -inf, -inf), Vec(inf, inf, inf));
	}
};

class Sphere : public Obj {
	public:
	Vec c;
	double r;

	Sphere(double r_= 0, Vec c_=0) { c=c_; r=r_; }
	double intersect(const Ray& ray) const {
		double b = ((ray.o-c)*2).dot(ray.d);
		double c_ = (ray.o-c).dot((ray.o-c)) - (r*r);
		double disc = b*b - 4*c_;
		if (disc<0) return 0;
		else disc = sqrt(disc);
		double sol1 = -b + disc;
		double sol2 = -b - disc;
		return (sol2>eps) ? sol2/2 : ((sol1>eps) ? sol1/2 : 0);
	}

	Vec normal(const Vec& p0) const {
		return (p0 - c).norm();
	}
	AABB getAABB() const {
		return AABB(Vec(c.x - r, c.y - r, c.z - r), Vec(c.x + r, c.y + r, c.z + r));
	}
};

class Intersection {
	public:
	Intersection() { t = inf; object = nullptr; }
	Intersection(double t_, Obj* object_) { t = t_; object = object_; }
	operator bool() { return object != nullptr; }
	double t;
	Obj* object;
};

class Scene {
	class BoundInfo {
	public:
		BoundInfo(Obj* object_)
		{
			object = object_;
			aabb = object->getAABB();
			centroid = (aabb.min + aabb.max) / 2.0;
		}

		Obj* object;
		AABB aabb;
		Vec centroid;
	};

	class BoundCentroidComparator {
	public:
		BoundCentroidComparator(size_t dimension_) { dimension = dimension_; }

		bool operator()(const BoundInfo& first, const BoundInfo& second) const {
			return first.centroid[dimension] < second.centroid[dimension];
		}

	private:
		size_t dimension;
	};

	class BoundingVolumeHierarchyNode {
	public:
		AABB aabb;
		union {
			uint32_t firstObjIndex;	// leaf node: index of the first Obj* (in boundedObjects)
			uint32_t secondChildIndex;	// interior node: index of the second child node (in boundingVolumeHierarchy)
		};
		uint8_t objectCount;	// leaf node: number of objects in boundedObjects (> 0); interior node: == 0
		uint8_t splitAxis;
	};

	public:
	void add(Obj* object) {
		AABB boundingBox = object->getAABB();
		if (boundingBox.unbounded()) {
			unboundedObjects.push_back(object);
		} else {
			// append the information needed for building the BVH to the boundInfos vector
			boundInfos.push_back(BoundInfo(object));
		}
	}

	void rebuildBVH(uint8_t maxObjectsPerLeaf) {
		if (maxObjectsPerLeaf == 0) {
			maxObjectsPerLeaf = 1;
		}

		// build the actual BVH, i.e. a depth-first representation stored in the boundingVolumeHierarchy vector
		boundingVolumeHierarchy.clear();
		boundingVolumeHierarchy.reserve(4 * boundInfos.size());	// this should be sufficient for well-balanced trees
		splitBoundsRecursively(boundInfos.begin(), boundInfos.end(), maxObjectsPerLeaf);

		// finally, grab just the pointers from boundInfos and place them in the boundedObjects vector
		//TODO: benchmark this to see if it's really faster than just reusing the boundInfos vector
		boundedObjects.clear();
		boundedObjects.reserve(boundInfos.size());
		for (auto& boundInfo : boundInfos) {
			boundedObjects.push_back(boundInfo.object);
		}
	}

	void splitBoundsRecursively(const std::vector<BoundInfo>::iterator begin, const std::vector<BoundInfo>::iterator end, uint8_t maxObjectsPerLeaf) {
		// append a node
		boundingVolumeHierarchy.push_back(BoundingVolumeHierarchyNode());
		std::vector<BoundingVolumeHierarchyNode>::iterator thisNode = boundingVolumeHierarchy.end(); --thisNode;
	
		const uint32_t objectCount = end - begin;
		if (objectCount <= maxObjectsPerLeaf) {
			// make thisNode a leaf
			//TODO: idea: we can still sort the entities along the largestDimension and store the splitAxis to enable front-to-back tracing later
			thisNode->objectCount = objectCount;
			for (auto iter = begin; iter != end; ++iter) {
				thisNode->aabb = enclose(thisNode->aabb, iter->aabb);
			}
			thisNode->firstObjIndex = begin - boundInfos.begin();
		} else {
			// make thisNode an interior node
			thisNode->objectCount = 0;
	
			AABB centroidBound;
			for (auto iter = begin; iter != end; ++iter) {
				centroidBound = enclose(centroidBound, iter->centroid);
			}
			thisNode->splitAxis = centroidBound.largestDimension();
			const auto medianIter = begin + (end - begin) / 2;
			std::nth_element(begin, medianIter, end, BoundCentroidComparator(thisNode->splitAxis));
	
			// recursively call for left and right child and note the index of the second child in-between
			const size_t firstChildIndex = boundingVolumeHierarchy.size();	// the first child follows this one directly (which is also why we don't store the index)
			splitBoundsRecursively(begin, medianIter, maxObjectsPerLeaf);
			thisNode->secondChildIndex = boundingVolumeHierarchy.size();	// the second child comes after the first child's tree
			splitBoundsRecursively(medianIter, end, maxObjectsPerLeaf);
	
			// the world bound of this node encloses the ones of both children
			thisNode->aabb = enclose(boundingVolumeHierarchy[firstChildIndex].aabb, boundingVolumeHierarchy[thisNode->secondChildIndex].aabb);
		}
	}

	Intersection intersect(const Ray& ray) const {
		Intersection closestIntersection;

		// first intersect all unbounded objects; these may reduce the intersection's t, which can speed up BVH traversal afterwards
		for (auto iter = unboundedObjects.begin(); iter != unboundedObjects.end(); ++iter) {
			double t = (*iter)->intersect(ray);
			if (t > eps && t < closestIntersection.t) {
				closestIntersection.t = t;
				closestIntersection.object = *iter;
			}
		}

		if (!boundingVolumeHierarchy.empty()) {
			Vec inverseDirection(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);	// precalculated for better performance

			uint32_t stack[64];	//TODO: make sure this depth is never exceeded while building the BVH
			uint32_t stackSize = 0;
			uint32_t nodeNumber = 0;
			while (true) {
				const BoundingVolumeHierarchyNode& node = boundingVolumeHierarchy[nodeNumber];
				if (node.aabb.intersect(ray, inverseDirection, closestIntersection.t)) {
					if (node.objectCount > 0) {	// leaf
						// intersect ray with the entities in the leaf
						for (size_t objectNumber = 0; objectNumber != node.objectCount; ++objectNumber) {
							double t = boundedObjects[node.firstObjIndex + objectNumber]->intersect(ray);
							if (t > eps && t < closestIntersection.t) {
								closestIntersection.t = t;
								closestIntersection.object = boundedObjects[node.firstObjIndex + objectNumber];
							}
						}
						if (stackSize == 0) {
							break;
						}
						nodeNumber = stack[--stackSize];
					} else {	// interior
						if (ray.d[node.splitAxis] < 0) {
							// stack the left node and try the right (closer) one next
							stack[stackSize++] = nodeNumber + 1;
							nodeNumber = node.secondChildIndex;
						} else {
							// stack the right node and try the left (closer) one next
							stack[stackSize++] = node.secondChildIndex;
							nodeNumber = nodeNumber + 1;
						}
					}
				} else {	// no intersection with this node; try next
					if (stackSize == 0) {
						break;
					}
					nodeNumber = stack[--stackSize];
				}
			}
		}

		return closestIntersection;
	}

	private:
	vector<BoundInfo> boundInfos;	// information needed during BVH build
	vector<BoundingVolumeHierarchyNode> boundingVolumeHierarchy;	// the nodes of the BVH tree (in depth-first order)
	vector<Obj*> boundedObjects;	// objects within the same BVH leaf node are placed next to each other
	vector<Obj*> unboundedObjects;	// e.g. planes; they don't go into the BVH
};

// Class for generating the Halton low-discrepancy series for Quasi
// Monte Carlo integration.
class Halton {
	double value, inv_base;
	public:
	void number(int i,int base) {
		double f = inv_base = 1.0 / base;
		value = 0.0;
		while(i > 0) {
			value += f * (double)(i%base);
			i /= base;
			f *= inv_base;
		}
	}
	void next() {
		double r = 1.0 - value - 0.0000001;
		if(inv_base<r) value += inv_base;
		else {
			double h = inv_base, hh;
			do {hh = h; h *= inv_base;} while(h >=r);
			value += hh + h - 1.0;
		}
	}
	double get() { return value; }
};

// Structures and classes needed for PSSMLT
class Path {
	struct PathPoint {
		float value; // u[i].value
		int modify; // u[i].modify
		stack<float> s; // history of values
	};
	vector<PathPoint> u;
	public:
	float getValue(int i) { return u[i].value; }
	int getModify(int i) { return u[i].modify; }
	int getSize() { return u.size(); }
	int getStackSize(int i) {return u[i].s.size();}
	void setValue(int i, float value) { u[i].value = value; }
	void setModify(int i, int modify) { u[i].modify = modify; }
	void add(float value, int modify) {
		PathPoint pathpt;
		pathpt.value = value;
		pathpt.modify = modify;
		stack<float> s;
		s.push(value);
		pathpt.s = s;
		u.push_back(pathpt);
	}
	Path() {
		add(RND2, 0);
		add(RND2, 0);
	}
	void saveState(int i, float value) {
		u[i].s.push(value);
	}
	void clearStack() {
		for (unsigned int i=0; i<u.size(); i++) {
			stack<float> s;
			u[i].s = s;
		}
	}
	void restoreState() {
		for (unsigned int i=0; i<u.size(); i++) {
			float value = u[i].value;
			while (!u[i].s.empty()) {
				value = u[i].s.top();
				u[i].s.pop();
			}
			u[i].value = value;
		}
	}
	void resetModify() {
		for (unsigned int i=0; i<u.size(); i++) {
			u[i].modify = 0;
		}
	}
};

struct Contrib {
	int x, y; // affected pixel
	Vec c; // contribution (RGB Vector)
};
struct Sample {
	float w; // weight of the sample
	Contrib c; // contribution of the sample
};

class MetroSampler {
	Path u;
	Sample oldsample;
	float oldI;
	float b;
	int M;
	int large_step;
	int time;
	int large_step_time;
	public:
	MetroSampler(const Path& u1, const float b1, const int M1, const int ls) 
	: u(u1), b(b1), M(M1), large_step(ls), time(ls), large_step_time(ls)
	{
		oldI = 0.0f;
		float w = 0;
		int x = 0; int y = 0;
		Vec c = Vec(0,0,0);
		oldsample.w = w;
		oldsample.c.x = x; oldsample.c.y = y;
		oldsample.c.c = c;
	}
	float Mutate( float value ) {
		float s1 = 1./1024, s2 = 1./64;
		float dv = s2*exp(-log(s2/s1)*RND2);
		if (RND2 < 0.5) {
			value += dv; if (value > 1) value -= 1;
		} else {
			value -= dv; if (value < 0) value += 1;
		}
		return value;
	}
	float PrimarySample(int i) {
		if (i>=u.getSize()) {
			// Then the path vector u is not big enough and we need to add random values
			// We add 4 at the same time to be sure (cause trace() uses at most 4 values)
			u.add(RND2, time);
			u.add(RND2, time);
			u.add(RND2, time);
			u.add(RND2, time);
		} else {
			if (u.getModify(i) < time) {
				if (large_step) { // large step
					u.saveState(i, u.getValue(i)); // save state
					u.setModify(i, time);
					u.setValue(i, RND2);
					if (time!=1) printf("\nDid a large step, i=%d\n",i);
				} else { // small step
					if (u.getModify(i) < large_step_time) {
						u.setModify(i, large_step_time);
						u.setValue(i, RND2);
						if (time!=1) printf("\nSimulated a large step, i=%d\n",i);
					}
					// lazy evaluation of mutations
					while (u.getModify(i) < time-1) {
						u.setValue(i, Mutate(u.getValue(i)));
						u.setModify(i, u.getModify(i)+1);
					}
					u.saveState(i, u.getValue(i)); // save state
					u.setValue(i, Mutate(u.getValue(i)));
					u.setModify(i, u.getModify(i)+1);
				}
			}
		}
		return u.getValue(i);
	}
	Sample Next(const float I, const Contrib& contrib) {
		float a = min(1.f, I/oldI); // accept prob.
		Sample newsample;
		Sample contribsample;
		newsample.c = contrib;
		newsample.w = (a+large_step)/(I/b+plarge)/M;
		// cumulate weight
		oldsample.w += (1-a)/(oldI/b+plarge)/M;
		if (RND2 < a) { // accept
			//printf("Accept, a=%1.3f\n",a);
			oldI = I;
			contribsample = oldsample;
			oldsample = newsample;
			if (large_step) large_step_time = time;
			time++;
			u.clearStack(); // no state restoration
		} else { // reject
			//printf("Reject, a=%1.3f\n",a);
			contribsample = newsample;
			u.restoreState();
		}
		large_step = (RND2 < plarge) ? 1 : 0;
		return contribsample;
	}
	void setPath(float value0, float value1) {
		// Method used to override evaluation of the 2 first values of the path
		// Only for seeds evaluation
		u.setValue(0, value0);
		u.setValue(1, value1);
		u.resetModify();
	}
	Path getPath() {
		// Returns current path
		// Only for seeds evaluation
		return u;
	}
};

float luminance(Vec clr){
	return (float)0.3*clr.x + 0.6*clr.y + 0.1*clr.z;
}

// Input is the pixel offset, output is the appropriate coordinate
// on the image plane
Vec camcr(const double x, const double y) {
	double w=width;
	double h=height;
	float fovx = PI/4;
	float fovy = (h/w) * fovx;
	return Vec(((2*x-w)/w) * tan(fovx),
				-((2*y-h)/h) * tan(fovy),
				-1.0);
}

// Uniform sampling on a hemisphere to produce outgoing ray directions.
// courtesy of http://www.rorydriscoll.com/2009/01/07/better-sampling/
Vec hemisphere(double u1, double u2) {
	const double r = sqrt(1.0-u1*u1);
	const double phi = 2 * PI * u2;
	return Vec(cos(phi)*r, sin(phi)*r, u1);
}

void trace(Ray &ray, const Scene& scene, int depth, Vec& clr, pl& params, MetroSampler& ms) {
	// Russian roulette: starting at depth 5, each recursive step will stop with a probability of 0.1
	double rrFactor = 1.0;
	if (depth >= 5) {
		const double rrStopProbability = 0.1;
		if (ms.PrimarySample(4*(depth+1)-2) <= rrStopProbability) {
			return;
		}
		rrFactor = 1.0 / (1.0 - rrStopProbability);
	}

	Intersection intersection = scene.intersect(ray);
	if (!intersection) return;

	// Travel the ray to the hit point where the closest object lies and compute the surface normal there.
	Vec hp = ray.o + ray.d * intersection.t;
	Vec N = intersection.object->normal(hp);
	ray.o = hp;
	// Add the emission, the L_e(x,w) part of the rendering equation, but scale it with the Russian Roulette
	// probability weight.
	const double emission = intersection.object->emission;
	clr = clr + Vec(emission, emission, emission) * rrFactor;

	// Diffuse BRDF - choose an outgoing direction with hemisphere sampling.
	if(intersection.object->type == 1) {
		Vec rotX, rotY;
		ons(N, rotX, rotY);
		Vec sampledDir = hemisphere(ms.PrimarySample(4*(depth+1)-1),ms.PrimarySample(4*(depth+1)+0));
		Vec rotatedDir;
		rotatedDir.x = Vec(rotX.x, rotY.x, N.x).dot(sampledDir);
		rotatedDir.y = Vec(rotX.y, rotY.y, N.y).dot(sampledDir);
		rotatedDir.z = Vec(rotX.z, rotY.z, N.z).dot(sampledDir);
		ray.d = rotatedDir;	// already normalized
		double cost=ray.d.dot(N);
		Vec tmp;
		trace(ray,scene,depth+1,tmp,params,ms);
		clr = clr + (tmp.mult(intersection.object->cl)) * cost * 0.1 * rrFactor;
	}

	// Specular BRDF - this is a singularity in the rendering equation that follows
	// delta distribution, therefore we handle this case explicitly - one incoming
	// direction -> one outgoing direction, that is, the perfect reflection direction.
	if(intersection.object->type == 2) {
		double cost = ray.d.dot(N);
		ray.d = (ray.d - N*(cost*2)).norm();
		Vec tmp = Vec(0,0,0);
		trace(ray,scene,depth+1,tmp,params,ms);
		clr = clr + tmp * rrFactor;
	}

	// Glass/refractive BRDF - we use the vector version of Snell's law and Fresnel's law
	// to compute the outgoing reflection and refraction directions and probability weights.
	if(intersection.object->type == 3) {
		double n = params["refr_index"];
		double R0 = (1.0-n)/(1.0+n);
		R0 = R0*R0;
		if(N.dot(ray.d)>0) { // we're inside the medium
			N = N*-1;
			n = 1/n;
		}
		n=1/n;
		double cost1 = (N.dot(ray.d))*-1; // cosine of theta_1
		double cost2 = 1.0 - n*n*(1.0-cost1*cost1); // cosine of theta_2
		double Rprob = R0 + (1.0-R0) * pow(1.0 - cost1, 5.0); // Schlick-approximation
		if (cost2 > 0 && ms.PrimarySample(4*(depth+1)+1) > Rprob) { // refraction direction
			ray.d = ((ray.d*n)+(N*(n*cost1-sqrt(cost2)))).norm();
		}
		else { // reflection direction
			ray.d = (ray.d+N*(cost1*2)).norm();
		}
		Vec tmp;
		trace(ray,scene,depth+1,tmp,params,ms);
		clr = clr + tmp * 1.15 * rrFactor;
	}
}

int main() {
	srand(time(NULL));
	pl params;
	Scene scene;
	auto add=[&scene](Obj* s, Vec cl, double emission, int type) {
			s->setMat(cl,emission,type);
			scene.add(s);
	};

	// Radius, position, color, emission, type (1=diff, 2=spec, 3=refr) for spheres
	add(new Sphere(1.05,Vec(-0.75,-1.45,-4.4)),Vec(4,8,4),0,2); // Middle sphere
//	add(new Sphere(0.45,Vec(0.8,-2.05,-3.7)),Vec(10,10,1),0,3); // Right sphere
	add(new Sphere(0.5,Vec(2.0,-2.05,-3.7)),Vec(10,10,1),0,3); // Right sphere
	add(new Sphere(0.6,Vec(-1.75,-1.95,-3.1)),Vec(4,4,12),0,1); // Left sphere
	// Position, normal, color, emission, type for planes
	add(new Plane(2.5,Vec(0,1,0)),Vec(6,6,6),0,1); // Bottom plane
	add(new Plane(5.5,Vec(0,0,1)),Vec(6,6,6),0,1); // Back plane
	add(new Plane(2.75,Vec(1,0,0)),Vec(10,2,2),0,1); // Left plane
	add(new Plane(2.75,Vec(-1,0,0)),Vec(2,10,2),0,1); // Right plane
	add(new Plane(3.0,Vec(0,-1,0)),Vec(6,6,6),0,1); // Ceiling plane
	add(new Plane(0.5,Vec(0,0,-1)),Vec(6,6,6),0,1); // Front plane
	add(new Sphere(0.5,Vec(0,1.9,-3)),Vec(0,0,0),5000,1); // Light
	scene.rebuildBVH(1);	// add all objects before this call

	params["refr_index"] = 1.5;
	params["spp"] = 4.0; // samples per pixel
	
	Vec **pix = new Vec*[width];
	for(int i=0;i<width;i++) {
		pix[i] = new Vec[height];
	}

	const double spp = params["spp"];
	// correlated Halton-sequence dimensions
	//Halton hal, hal2;
	//hal.number(0,2);
	//hal2.number(0,2);

	clock_t start = clock();
	
	// Initialisation phase for PSSMLT
	printf("\nGenerating path seeds...\n");
	// We can instantiate a MetroSampler with dummy values, since we won't use the Next() function
	// and by default we make only large steps
	Path u;
	MetroSampler msi(u, 1, 1, 1);
	// Initialization of the variables to compute b and u1
	float b = 0.0f;
	Path u1;
	float I = 0;
	float bestI = 0;
	// We create as many seed paths as we have pixels in the image
	for (int col=0;col<width;col++) {
		for(int row=0;row<height;row++) {
			Vec color;
			Ray ray;
			ray.o = (Vec(0,0,0)); // rays start out from here
			Vec cam = camcr(col,row); // construct image plane coordinates
			// sets the 2 first coordinates of the seed path
			msi.setPath((float)col/(float)width, (float)row/(float)height);
			ray.d = (cam - ray.o).norm(); // point from the origin to the camera plane
			trace(ray,scene,0,color,params,msi);
			I = luminance(color);
			b += I;
			if (I > bestI) {
				bestI = I;
				u1 = msi.getPath();
			}
		}
	}
	printf("Init OK\n");
	b /= height*width;
	u1.resetModify();
	int M = (int)spp;
	
	// Now we have the unbiaised initialization values, we can instantiate the real MetroSampler
	MetroSampler ms(u1, b, M, 0);

	// Main loop
	//#pragma omp parallel for schedule(dynamic)
	for (int col=0;col<width;col++) {
		fprintf(stdout,"\rRendering: %1.0fspp %8.2f%%",spp,(double)col/width*100);
		for(int row=0;row<height;row++) {
			for(int s=0;s<spp;s++) {
				//fprintf(stdout,"Rendering: %1.0fspp %8.2f%% \n",spp,(double)col/width*100);
				Vec color;
				Ray ray;
				ray.o = (Vec(0,0,0)); // rays start out from here
				// construct image plane coordinates from path first 2 values of sample
				double x = (double)(ms.PrimarySample(0)*(width-1));
				double y = (double)(ms.PrimarySample(1)*(height-1));
				Vec cam = camcr(x, y);
				ray.d = (cam - ray.o).norm(); // point from the origin to the camera plane
				trace(ray,scene,0,color,params,ms);
				// Computing contribution
				I = luminance(color);
				Contrib contribution;
				contribution.x = (int)round(x);
				contribution.y = (int)round(y);
				contribution.c = color;
				//printf("x=%d  y=%d  I=%2.2f\n", contribution.x, contribution.y, I);
				// Call the Metropolis Sampler
				Sample sample = ms.Next(I, contribution);
				// Compute contribution of the sample
				double w = (double)sample.w;
				int j1 = sample.c.x; int j2 = sample.c.y;
				Vec clr = sample.c.c;
				//I = luminance(clr);
				//printf("j1=%d  j2=%d  w=%2.3f I=%2.2f\n", j1, j2, w, I);
				pix[j1][j2] = pix[j1][j2] + clr * w; // write the contributions
			}
		}
	}

	FILE *f = fopen("ray_PSSMLT.ppm", "w");
	fprintf(f, "P3\n%d %d\n%d\n ",width,height,255);
	for (int row=0;row<height;row++) {
		for (int col=0;col<width;col++) {
			fprintf(f,"%d %d %d ", min((int)pix[col][row].x,255), min((int)pix[col][row].y,255), min((int)pix[col][row].z,255));
		}
		fprintf(f, "\n");
	}
	fclose(f);
	clock_t end = clock();
	double t = (double)(end-start)/CLOCKS_PER_SEC;
	printf("\nRender time: %fs.\n",t);
	return 0;
}
