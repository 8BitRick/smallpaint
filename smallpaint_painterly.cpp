// smallpaint by karoly zsolnai - zsolnai@cg.tuwien.ac.at
//
// compilation by: g++ smallpaint.cpp -O3 -std=gnu++0x -fopenmp
// uses >=gcc-4.5.0
// render, modify, create new scenes, tinker around, and most of all:
// have fun!
//
// If you have problems/DLL hell with the compilation under Windows and
// MinGW, you can add these flags: -static-libgcc -static-libstdc++
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

#if defined(__linux) || defined(__linux__)
	unsigned int seed = time(NULL);
	#define RND (2.0*(double)rand_r(&seed)/RAND_MAX-1.0)
	#define RND2 ((double)rand_r(&seed)/RAND_MAX)
#endif

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
	#define RND (2.0*(double)rand()/RAND_MAX-1.0)
	#define RND2 ((double)rand()/RAND_MAX)
#endif

#define PI 3.1415926536
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const int width=900, height=900;
const double inf=1e9;
const double eps=1e-4;
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
	Vec operator%(Vec &b){return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
};

struct Ray {
	Vec o, d;
	Ray(Vec o0 = 0, Vec d0 = 0) { o = o0, d = d0.norm(); }
};

class Obj {
	public:
	Vec cl;
	double emission;
	int type;
	void setMat(Vec cl_=0, double emission_=0, int type_=0) { cl=cl_; emission=emission_; type=type_; }
	virtual double intersect(const Ray&) const =0;
	virtual Vec normal(const Vec&) const =0;
};

class Plane : public Obj {
	public:
	Vec n;
	double d;
	Plane(double d_ =0, Vec n_=0) { d=d_; n=n_; }
	double intersect(const Ray& ray) const {
		double d0 = n.dot(ray.d);
		if(d0!=0) {
			double t = -1*(((n.dot(ray.o))+d)/d0);
			return (t>eps) ? t : 0;
		}
		else return 0;
	}
	Vec normal(const Vec& p0) const { return n; }
};

class Sphere : public Obj {
	public:
	Vec c;
	double r;

	Sphere(double r_= 0, Vec c_=0) { c=c_; r=r_; }
	double intersect(const Ray& ray) const {
		double b = ((ray.o-c)*2).dot(ray.d);
		double c_ = (ray.o-c).dot((ray.o-c))-(r*r);
		double disc = b*b - 4*c_;
		if (disc<0) return 0;
		else disc = sqrt(disc);
		double sol1 = -b+disc;
		double sol2 = -b-disc;
		return (sol2>eps) ? sol2/2 : ((sol1>eps) ? sol1/2 : 0);
	}

	Vec normal(const Vec& p0) const {
		return Vec((p0.x-c.x)/r,
					(p0.y-c.y)/r,
					(p0.z-c.z)/r);
	}
};

class Halton {
	double value, inv_base;
	public:
	void number(int i,int base) {
		double f = inv_base = 1.0/base;
		value = 0.0;
		while(i>0) {
			value += f*(double)(i%base);
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

Vec camcr(const double x, const double y) {
	double w=width;
	double h=height;
	float fovx = PI/4;
	float fovy = (h/w)*fovx;
	return Vec(((2*x-w)/w)*tan(fovx),
				((2*y-h)/h)*tan(fovy),
				-1.0);
}

// a messed up sampling function (at least in this context).
// courtesy of http://www.rorydriscoll.com/2009/01/07/better-sampling/
Vec hemisphere(double u1, double u2) {
	const double r=sqrt(1.0-u1*u1);
	const double phi=2*PI*u2;
	return Vec(cos(phi)*r,sin(phi)*r,u1);
}

void trace(Ray &ray, const vector<Obj*>& scene, int depth, Vec& clr, pl& params, Halton& hal, Halton& hal2) {
	if (depth >= 20) return;
	double t;
	int id = -1;
	double mint = inf;
	const unsigned int size = scene.size();

			for(unsigned int x=0;x<size;x++) {
					t = scene[x]->intersect(ray);
					if (t>eps && t < mint) { mint=t; id=x; }
				}
			if (id == -1) return;
			Vec hp = ray.o + ray.d*mint;
			Vec N = scene[id]->normal(hp);
			ray.o = hp;
			clr = clr + Vec(scene[id]->emission,scene[id]->emission,scene[id]->emission)*2;

			if(scene[id]->type == 1) {
				hal.next();
				hal2.next();
				ray.d = (N+hemisphere(hal.get(),hal2.get()));
				double cost=ray.d.dot(N);
				Vec tmp=Vec();
				trace(ray,scene,depth+1,tmp,params,hal,hal2);
				clr.x += cost*(tmp.x*scene[id]->cl.x)*0.1;
				clr.y += cost*(tmp.y*scene[id]->cl.y)*0.1;
				clr.z += cost*(tmp.z*scene[id]->cl.z)*0.1;
			}

			if(scene[id]->type == 2) {
				double cost = ray.d.dot(N);
				ray.d = (ray.d-N*(cost*2)).norm();
				Vec tmp = Vec(0,0,0);
				trace(ray,scene,depth+1,tmp,params,hal,hal2);
				clr = clr + tmp;
			}

			if(scene[id]->type == 3) {
				double n = params["refr_index"];
				if(N.dot(ray.d)>0) {
					N = N*-1;
					n = 1/n;
				}
				n=1/n;
				double cost1 = (N.dot(ray.d))*-1;
				double cost2 = 1.0-n*n*(1.0-cost1*cost1);
				if(cost2 > 0) {
					ray.d = (ray.d*n)+(N*(n*cost1-sqrt(cost2)));
					ray.d = ray.d.norm();
					Vec tmp = Vec(0,0,0);
					trace(ray,scene,depth+1,tmp,params,hal,hal2);
					clr = clr + tmp;
				}
				else return;
			}
}
int main() {

	srand(time(NULL));
	pl params;
	vector<Obj*> scene;
	auto add=[&scene](Obj* s, Vec cl, double emission, int type) {
			s->setMat(cl,emission,type);
			scene.push_back(s);
	};

	// Radius, position, color, emission, type (1=diff, 2=spec, 3=refr) for spheres
	add(new Sphere(1.05,Vec(1.45,-0.75,-4.4)),Vec(4,8,4),0,2); // Middle sphere
	add(new Sphere(0.45,Vec(2.05,0.8,-3.7)),Vec(10,10,1),0,3); // Right sphere
	add(new Sphere(0.6,Vec(1.95,-1.75,-3.1)),Vec(4,4,12),0,1); // Left sphere
	// Position, normal, color, emission, type for planes
	add(new Plane(2.5,Vec(-1,0,0)),Vec(6,6,6),0,1); // Bottom plane
	add(new Plane(5.5,Vec(0,0,1)),Vec(6,6,6),0,1); // Back plane
	add(new Plane(2.75,Vec(0,1,0)),Vec(10,2,2),0,1); // Left plane
	add(new Plane(2.75,Vec(0,-1,0)),Vec(2,10,2),0,1); // Right plane
	add(new Plane(3.0,Vec(1,0,0)),Vec(6,6,6),0,1); // Ceiling plane
	add(new Plane(0.5,Vec(0,0,-1)),Vec(6,6,6),0,1); // Front plane
	const Vec L0 = (Vec(-1.9,0,-3));
	add(new Sphere(eps,L0),Vec(0,0,0),120,1); // Light
	params["refr_index"] = 1.9;
	params["spp"] = 8.0; // samples per pixel

	Vec **pix = new Vec*[width];
	for(int i=0;i<width;i++) {
			pix[i] = new Vec[height];
	}

	const double spp = params["spp"];
	// correlated Halton-sequence dimensions
	Halton hal, hal2;
	hal.number(0,2);
	hal2.number(0,2);

	clock_t start = clock();

	#pragma omp parallel for schedule(dynamic) firstprivate(hal,hal2)
	for (int i=0;i<width;i++) {
		fprintf(stdout,"\rRendering: %1.0fspp %8.2f%%",spp,(double)i/width*100);
		for(int j=0;j<height;j++) {
			for(int s=0;s<spp;s++) {
				Vec c;
				Ray ray;
				ray.o = (Vec(0,0,0));
				Vec cam = camcr(i,j);
				cam.x = cam.x + RND/700;
				cam.y = cam.y + RND/700;
				ray.d = (cam - ray.o).norm();
				trace(ray,scene,0,c,params,hal,hal2);
				pix[i][j].x += c.x*(1/spp);
				pix[i][j].y += c.y*(1/spp);
				pix[i][j].z += c.z*(1/spp);
			}
		}
	}

	FILE *f = fopen("ray.ppm", "w");
	fprintf(f, "P3\n%d %d\n%d\n ",width,height,255);
	for (int i=0;i<width;i++)
	for (int j=0;j<height;j++) {
		fprintf(f,"%d %d %d ", min((int)pix[i][j].x,255), min((int)pix[i][j].y,255), min((int)pix[i][j].z,255));
	}
	fclose(f);
	clock_t end = clock();
	double t = (double)(end-start)/CLOCKS_PER_SEC;
	printf("\nRender time: %fs.\n",t);
	return 0;
}
