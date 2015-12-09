#ifndef Ray_hpp
#define Ray_hpp

#include "Vec.hpp"

// Rays have origin and direction.
// The direction vector should always be normalized.
struct Ray {
	Vec o, d;
	inline Ray(Vec o0 = 0, Vec d0 = 0) { o = o0, d = d0.norm(); }
};

#endif
