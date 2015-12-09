#ifndef AABB_hpp
#define AABB_hpp

#include <algorithm>

#include "Ray.hpp"
#include "constants.hpp"

struct AABB {
	inline AABB() { min = Vec(inf, inf, inf); max = Vec(-inf, -inf, -inf); }	// an empty interval
	inline AABB(Vec min_, Vec max_) { min = min_; max = max_; }
	inline bool unbounded() const { return min.x == -inf || min.y == -inf || min.z == -inf || max.x == inf || max.y == inf || max.z == inf; }
	inline size_t largestDimension() const {
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
	inline bool intersect(const Ray& ray, const Vec& inverseDirection, double closestKnownT) const {
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

inline AABB enclose(const AABB& firstBoundingBox, const AABB& secondBoundingBox) {
	AABB ret;

	ret.min.x = std::min(firstBoundingBox.min.x, secondBoundingBox.min.x);
	ret.min.y = std::min(firstBoundingBox.min.y, secondBoundingBox.min.y);
	ret.min.z = std::min(firstBoundingBox.min.z, secondBoundingBox.min.z);

	ret.max.x = std::max(firstBoundingBox.max.x, secondBoundingBox.max.x);
	ret.max.y = std::max(firstBoundingBox.max.y, secondBoundingBox.max.y);
	ret.max.z = std::max(firstBoundingBox.max.z, secondBoundingBox.max.z);

	return ret;
}

inline AABB enclose(const AABB& boundingBox, const Vec& point) {
	AABB ret;

	ret.min.x = std::min(boundingBox.min.x, point.x);
	ret.min.y = std::min(boundingBox.min.y, point.y);
	ret.min.z = std::min(boundingBox.min.z, point.z);

	ret.max.x = std::max(boundingBox.max.x, point.x);
	ret.max.y = std::max(boundingBox.max.y, point.y);
	ret.max.z = std::max(boundingBox.max.z, point.z);

	return ret;
}

#endif
