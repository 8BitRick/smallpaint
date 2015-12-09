// smallpaint by karoly zsolnai - zsolnai@cg.tuwien.ac.at
// extended for heterogeneous participating media by georg sperl - e1025854@student.tuwien.ac.at
//
// compilation by: g++ smallmedia.cpp -O3 -std=gnu++0x -fopenmp
// uses >=gcc-4.8.1
//
// For more informations on equiangular sampling, see:
// https://www.solidangle.com/research/egsr2012_volume.pdf
// https://www.shadertoy.com/view/Xdf3zB
// http://cg.tuwien.ac.at/~zsolnai/gfx/volumetric-path-tracing-with-equiangular-sampling-in-2k/
//
// Possible extensions:
// i) Multiple importance sampling of equiangular sampling and delta tracking as a sampling method
// ii) Residual ratio tracking

#include "smallmedia.hpp"

typedef unordered_map<string, double> pl; // parameter list

// Recursive tracing method
Vec trace(Ray &ray, const Scene& scene, int depth, pl& params) {

	// Maximum depth
	if (params["max_depth"] >= 0 && depth > params["max_depth"]) return Vec(0.0);

	// Russian Roulette
	double rrFactor = 1.0;
	if (params["rr_depth"] >= 0 && depth >= params["rr_depth"]) { //if russian roulette is activated and the starting depth is reached
		const double rrStopProbability = params["rr_prob"];
		if (RND2 <= rrStopProbability) {
			return Vec();
		}
		rrFactor = 1.0 / (1.0 - rrStopProbability);
	}

	// Traced color
	Vec clr;

	// Compute nearest surface intersection
	Intersection ints;
#ifdef SURFACES
	ints = scene.intersect(ray);
#endif

	double t_surface = ints.t; //inf if nothing was hit

	// Compute scattering event position
	double u = RND2;
	double t_scatter, pdf; //scattering event position on ray and probability
	//sampleUniform(u, inf, t_scatter, pdf); // Uniform sampling
	//sampleExponential(u, inf, scene.medium->sigma_t, t_scatter, pdf); // Exponential sampling
	sampleEquiAngular(u, inf, ray, scene.light->samplePos(), t_scatter, pdf); // Equiangular sampling

	Vec sp = ray.o + ray.d * t_scatter;

	// Sample the medium if the scattering event is nearer than the surface and the scattering position really is inside the medium
	if (t_scatter > eps && t_scatter < t_surface && scene.medium->inbounds(sp))
	{
		Vec l_pos = scene.light->samplePos(); // Light position

		Vec inscattering;

		// Calculate scattering coefficient at the scattering event position
		Vec sigma_s = scene.medium->sigma_s(sp);
		if (sigma_s.x > 0 && sigma_s.y > 0 && sigma_s.z > 0) // if the scattering coefficient is 0, inscattering does not contribute at all
		{
			// Explicit light sampling
			Vec L = l_pos - sp; // direction from scattering event to light
			double d = L.length(); // distance to light
			L.norm();

			Ray l_ray(sp, L);

			Intersection l_ints;
#ifdef SURFACES
			// cast shadow ray
			l_ints = scene.intersect(l_ray);
#endif
			double l_t = l_ints.t; // Distance to nearest intersection in light direction


			if (l_t >= d) // not shadowed (nearest intersection farther away than light)
			{
				double geom = 1.0 / (d * d); // quadr. attenuation
				//transmittance from scattering event to light
#ifdef RM
				Vec transmit = scene.medium->transmittanceRM(l_ray,0.0,d); // ray marching
#else
				Vec transmit = scene.medium->transmittanceDelta(l_ray,0,d, params["delta_samples"]); // delta tracking
#endif
#ifdef HOM
				transmit = scene.medium->transmittanceHom(d); // exponential
#endif
				inscattering += scene.light->emission() * INV_4PI * geom * transmit; //Note: /4pi because explicit light sampling chooses the one direction towards the light source
			}



#ifdef MULTIPLE_SCATTERING
				// trace a scatter ray with random direction and add the color to inscattering term
				Ray scatterRay(sp, sampleSphere(RND2, RND2));
				inscattering += trace(scatterRay, scene, depth + 1, params);
#endif
		}

		Vec emission = scene.medium->emission(sp) * INV_4PI * t_scatter; //emission of the media from ray.o to scattering event

		//transmittance from ray.o to the scattering event
#ifdef RM
		Vec transmit = scene.medium->transmittanceRM(ray, 0.0, t_scatter); // ray marching
#else
		Vec transmit = scene.medium->transmittanceDelta(ray, 0, t_scatter, params["delta_samples"]); // delta tracking
#endif
#ifdef HOM
		transmit = scene.medium->transmittanceHom(t_scatter); // exponential
#endif
		

		clr += (sigma_s * inscattering + emission) * transmit / pdf; // radiative transfer equation

	}
	else if (t_surface > eps && t_surface < inf) // medium not sampled -> sample surface
	{

		Vec hp = ray.o + ray.d * t_surface; // surface hit point
		Vec N = ints.object->normal(hp); // normal at the hit point

		Vec surface_clr;

		// explicit light sampling (only for diffuse surfaces)
		if (ints.object->type == 1)
		{
			Vec l_pos = scene.light->samplePos();
			Vec L = l_pos - hp;
			double d = L.length();
			L.norm();
			Ray l_ray(hp, L);
			Intersection l_ints = scene.intersect(l_ray); // cast a shadow ray towards light source

			if (l_ints.t > eps && l_ints.t < inf) //if the shadow ray "accidentally" hits an emissive surface the surface contributes its emission
			{
				const double emission = l_ints.object->emission;
				surface_clr += Vec(emission, emission, emission);
			}

			if (l_ints.t >= d) // not shadowed diffuse surface
			{
				//transmittance from surface to light
#ifdef RM
				Vec transmit = scene.medium->transmittanceRM(l_ray,0.0,d); // ray marching
#else
				Vec transmit = scene.medium->transmittanceDelta(l_ray, 0, d, params["delta_samples"]); // delta tracking
#endif
#ifdef HOM
				transmit = scene.medium->transmittanceHom(d); // exponential
#endif

				double nDotL = N.dot(L);
				double geom = nDotL / (d * d); //geometric term, including lambertian factor

				Vec expl_light = scene.light->emission() * INV_4PI * geom * transmit;
				expl_light = expl_light * ints.object->clr; //diffuse surface -> colour
				surface_clr += expl_light;
			}
		}

		// Surface interaction (reflection/refraction)

		// Diffuse BRDF - choose an outgoing direction with hemisphere sampling.
		if (ints.object->type == 1) {
			Vec rotX, rotY;
			ons(N, rotX, rotY); //construct a orthonormal system around the normal
			Vec sampledDir = sampleHemisphere(RND2, RND2);
			Vec rotatedDir;
			rotatedDir.x = Vec(rotX.x, rotY.x, N.x).dot(sampledDir); //transform the (sampled) hemisphere with the constructed ons
			rotatedDir.y = Vec(rotX.y, rotY.y, N.y).dot(sampledDir);
			rotatedDir.z = Vec(rotX.z, rotY.z, N.z).dot(sampledDir);
			Ray next(hp, rotatedDir);
			double cost = rotatedDir.dot(N); // attenuation (lambert) of reflected ray
			Vec tmp = trace(next, scene, depth + 1, params); // trace reflected ray
			surface_clr += tmp * ints.object->clr * cost;
		}

		// Specular BRDF - this is a singularity in the rendering equation that follows
		// delta distribution, therefore we handle this case explicitly - one incoming
		// direction -> one outgoing direction, that is, the perfect reflection direction.
		if (ints.object->type == 2) {
			double cost = ray.d.dot(N);
			Ray next;
			next.o = hp;
			next.d = (ray.d - N*(cost * 2)).norm();
			Vec tmp = trace(next, scene, depth + 1, params);
			surface_clr += tmp;
		}

		// Glass/refractive BRDF - we use the vector version of Snell's law and Fresnel's law
		// to compute the outgoing reflection and refraction directions and probability weights.
		if (ints.object->type == 3) {
			double n = params["refr_index"];
			double R0 = (1.0 - n) / (1.0 + n);
			R0 = R0*R0;
			if (N.dot(ray.d)>0) { // we're inside the refractive medium -> invert normal and swap n1/n2 (n1 = 1 = air, n2 = n) with n2/n1
				N = N*-1;
				n = 1 / n;
			}
			n = 1 / n; // used in cost2 as n*n to be 1/n^2
			double cost1 = (N.dot(ray.d))*-1; // cosine of theta1
			double cost2 = 1.0 - n*n*(1.0 - cost1*cost1); // cosine of theta_2 (squared)
			double Rprob = R0 + (1.0 - R0) * pow(1.0 - cost1, 5.0); // Schlick-approximation (probability of refraction)
			Ray next;
			next.o = hp;
			if (cost2 > 0 && RND2 > Rprob) { // refraction direction
				next.d = ((ray.d*n) + (N*(n*cost1 - sqrt(cost2)))).norm();
			}
			else { // reflection direction
				next.d = (ray.d + N*(cost1 * 2)).norm();
			}
			Vec tmp = trace(next, scene, depth + 1, params);
			surface_clr += tmp;
		}

		//transmittance from ray.o to surface
#ifdef RM
		Vec transmit = scene.medium->transmittanceRM(ray, 0.0, t_surface); // ray marching
#else
		Vec transmit = scene.medium->transmittanceDelta(ray, 0, t_surface, params["delta_samples"]); // delta tracking
#endif
#ifdef HOM
		transmit = scene.medium->transmittanceHom(t_surface); // exponential
#endif
		// Add the transmitted surface color
		clr += surface_clr * transmit;
	}

	// return accumulated color weighted with the russian roulette weight
	return clr * rrFactor;
}


int main() {
	srand(time(NULL));

	// SCENE SETUP:

	// Light (pos, intensity, color)
#ifndef SURFACES
	PointLight light(Vec(0.0, 0.0, -3.0), 15, Vec(1, 1, 1)); // Light in the center
#else
	PointLight light(Vec(0, 1.0, -3), 40, Vec(1, 1, 1));
#endif

	// Medium
	Vec sigma_a = Vec(0.1);				// absorption coefficient
	Vec sigma_s = Vec(0.05,0.10,0.07);	// scattering coefficient
	Vec emission = Vec(0);				// volumetric emission

	Vec position = Vec(-3, -3, -15);
	Vec size = Vec(6, 6, 30);

	// Density grid data
	const unsigned int n = 25; //Minimum 2 // less than 50..
	double densities[n*n*n];
	int type = 3; //0..uniform random, 1..checkerboard, 2..gradient, 3..perlin, 4..constant
	Grid grid(n, densities);
	fillGrid(type, grid);

	Medium medium(sigma_a, sigma_s, emission, grid, position, size);

	Scene scene;
	scene.light = &light;
	scene.medium = &medium;

	auto add = [&scene](Obj* s, Vec cl, double emission, int type) {
		s->setMat(cl, emission, type);
		scene.add(s);
	};

#ifdef CORNELLBOX
	// Cornell-Box
	add(new Plane(2.5, Vec(0, 1, 0)), Vec(1), 0, 1); // Bottom plane
	add(new Plane(5.5, Vec(0, 0, 1)), Vec(1), 0, 1); // Back plane
	add(new Plane(2.75, Vec(1, 0, 0)), Vec(1, .2, .2), 0, 1); // Left plane
	add(new Plane(2.75, Vec(-1, 0, 0)), Vec(.2, 1, .2), 0, 1); // Right plane
	add(new Plane(3.0, Vec(0, -1, 0)), Vec(1), 0, 1); // Ceiling plane
	add(new Plane(0.5, Vec(0, 0, -1)), Vec(1), 0, 1); // Front plane
#endif

	// Other objects
	// add( new Sphere(Radius, position, color, emission, type (1=diff, 2=spec, 3=refr)) for spheres
	add(new Sphere(1.05, Vec(-0.75, -1.45, -4.4)), Vec(0, 0, 0), 0, 2); // Middle sphere
	add(new Sphere(0.5, Vec(1.2, 0.6, -2.7)), Vec(0, 0, 0), 0, 3); // Right sphere
	add(new Sphere(0.6, Vec(-1.75, 1.2, -3.1)), Vec(0.3, 0.3, 1), 0, 1); // Left sphere

	// TRACE PARAMETERS:
	pl params;
	params["spp"] = 8; // samples per pixel
	params["refr_index"] = 2.3; // glass

	params["max_depth"] = 5;	// maximum trace depth (-1 = no maximum depth --> need to use RR!)
	params["rr_depth"] = 4;		// starting depth of russian roulette termination (-1 = no RR)
	params["rr_prob"] = 0.1;	// RR stopping probability

	params["delta_samples"] = 10; // delta tracking samples

	// RENDER PARAMETERS:
	bool gammacorrection = false;

	// Output image dimensions
	const int width = 480, height = 480;
	// render stencil (only renders the selected portion of the image)
	double minx = 0.0 * width;
	double maxx = 1.0 * width;
	double miny = 0.0 * height;
	double maxy = 1.0 * height;

	// RENDER SCENE //
	Vec **pix = new Vec*[width];
	for (int i = 0; i<width; i++) {
		pix[i] = new Vec[height];
	}

	const double spp = params["spp"];

	clock_t start = clock();

	#pragma omp parallel for schedule(dynamic)
	for (int col = 0; col<width; col++) {
		fprintf(stdout, "\rRendering: %1.0fspp %8.2f%%", spp, (double)col / width * 100);
		for (int row = 0; row<height; row++) {
			for (int s = 0; s<spp; s++) {
				Ray ray;
				ray.o = (Vec(0, 0, 0)); // rays start out from here
				Vec cam = camcr(col, row, width, height); // construct image plane coordinates
				cam.x = cam.x + RND / 700; // anti-aliasing for free (different rays per pixel)
				cam.y = cam.y + RND / 700;
				ray.d = (cam - ray.o).norm(); // point from the origin to the camera plane
				if (col >= minx && col < maxx && row >= miny && row < maxy) // render stencil area
				{
					Vec color = trace(ray, scene, 0, params);
					pix[col][row] = pix[col][row] + color / spp; // write the contributions
				}
			}
			//assert(pix[col][row].min() >= 0, "Negative final color: " + std::to_string(pix[col][row].min()));
		}
	}

	// SAVE IMAGE FILE //
	FILE *f = fopen("ray.ppm", "w");
	fprintf(f, "P3\n%d %d\n%d\n ", width, height, 255);
	for (int row = 0; row<height; row++) {
		for (int col = 0; col<width; col++) {
			fprintf(f, "%d %d %d ", clrToInt(pix[col][row].x, gammacorrection), clrToInt(pix[col][row].y, gammacorrection), clrToInt(pix[col][row].z, gammacorrection));
		}
		fprintf(f, "\n");
	}
	fclose(f);
	clock_t end = clock();
	double t = (double)(end - start) / CLOCKS_PER_SEC;
	printf("\nRender time: %fs.\n", t);
	return 0;
}
