// Multiple importance sampling demonstration in 100 lines of code (plus printing the results).
// It was made for the Rendering course at TU Wien:
// http://cg.tuwien.ac.at/courses/Rendering/
// Compiling this requires c++0x/c++11!
// Károly Zsolnai - zsolnai@cg.tuwien.ac.at
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <random>

#if defined(__linux) || defined(__linux__)
	unsigned int seed = time(NULL);
	#define RND (2.0*(double)rand_r(&seed)/RAND_MAX-1.0)
	#define RND2 ((double)rand_r(&seed)/RAND_MAX)
#endif

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
	#define RND (2.0*(double)rand()/RAND_MAX-1.0)
	#define RND2 ((double)rand()/RAND_MAX)
#endif

#define PI 3.1415926535

double f(double x) {
	double mu = PI/2;
	double sigma = PI/12;
	return ( 1 / ( sigma * sqrt(2*PI) ) ) * expf(- ( ((x-mu)*(x-mu)) / (2*sigma*sigma) ));
}

double gauss_pdf(double x, double mu, double sigma) {
	// Probability density function of the 1D normal distribution 
	// - http://en.wikipedia.org/wiki/Normal_distribution
	return ( 1 / ( sigma * sqrt(2*PI) ) ) * expf(- ( ((x-mu)*(x-mu)) / (2*sigma*sigma) ));
}

int main() {
	srand(time(NULL));
	printf("\nIntegrating f(x) = 1/(sigma*sqrt(2pi)) * e^(-(x-mu)^2/2sigma^2) from 0 to pi...\n\n");
	double sum1 = 0;
	double sum2 = 0;
	double MISsum1 = 0, MISsum2 = 0, MISsum3 = 0, MISsum4 = 0;
	int samples = 1000001;
	int report = 0;
	double result = 1; // the result of the integration to check the error against.
	double a = 0; // lower bound of integration
	double b = PI; // upper bound of integration
	double mu = (b-a)/2; // mean (expected value) of the normal distribution
	double sigma = (b-a)/12; // standard deviation of the normal distribution
	double f_x1, f_x2; // integrand
	double p_x1 = 1/(b-a); // Uniform distribution pdf
	double w1, w2; // MIS weights
	std::normal_distribution<> normal(mu, sigma); // c++11 normal distribution generator
	std::mt19937 mt(42);

	for(int i=0;i<samples;i++) {
		f_x1 = f(RND2*(b-a)); // sample the integrand with uniform distribution between the bounds
		double xi;
		// here we sample the integrand and make sure the sample is in the integration bounds
		for(;;) { xi = normal(mt); if(xi > a && xi < b) break; }
		f_x2 = f(xi);
		double p_x2 = gauss_pdf(xi, mu, sigma);
		// build the regular Monte Carlo estimators for both techniques
		sum1 += f_x1 / p_x1;
		sum2 += f_x2 / p_x2;

		// Eric Veach's excellent thesis contains lots of information on what is happening here:
		// http://graphics.stanford.edu/papers/veach_thesis/
 		// build the MIS estimators to combine the two techniques together
		w1 = 0.5; // naive averaging
		w2 = 0.5;
		MISsum1 += w1 * (f_x1 / p_x1) + w2 * (f_x2 / p_x2);

		w1 = p_x1 / (p_x1 + p_x2); // balance heuristic
		w2 = p_x2 / (p_x1 + p_x2);
		MISsum2 += w1 * (f_x1 / p_x1) + w2 * (f_x2 / p_x2);

		int beta = 2;
		w1 = pow(p_x1,beta) / (pow(p_x1,beta) + pow(p_x2,beta)); // power heuristic
		w2 = pow(p_x2,beta) / (pow(p_x1,beta) + pow(p_x2,beta)); 
		MISsum3 += w1 * (f_x1 / p_x1) + w2 * (f_x2 / p_x2);

		int beta_2 = 500;
		w1 = pow(p_x1,beta_2) / (pow(p_x1,beta_2) + pow(p_x2,beta_2)); // power heuristic
		w2 = pow(p_x2,beta_2) / (pow(p_x1,beta_2) + pow(p_x2,beta_2)); 
		MISsum4 += w1 * (f_x1 / p_x1) + w2 * (f_x2 / p_x2);

		if((int)log10(i) == report) {
			double estimator1 = sum1/i;
			double estimator2 = sum2/i;
			double MISestimator1 = MISsum1/i;
			double MISestimator2 = MISsum2/i;
			double MISestimator3 = MISsum3/i;
			double MISestimator4 = MISsum4/i;
			double error_estimator1 = fabs(result-estimator1);
			double error_estimator2 = fabs(result-estimator2);
			double error_MISestimator1 = fabs(result-MISestimator1);
			double error_MISestimator2 = fabs(result-MISestimator2);
			double error_MISestimator3 = fabs(result-MISestimator3);
			double error_MISestimator4 = fabs(result-MISestimator4);
			printf("#1 Monte Carlo estimator after %d samples is : %f, error: %f\n",i,estimator1,error_estimator1);
			printf("#2 Monte Carlo estimator after %d samples is : %f, error: %f\n",i,estimator2,error_estimator2);
			printf("#3 MIS combined estimator after %d samples is : %f, error: %f (Naive Averaging)\n",i,MISestimator1,error_MISestimator1);
			printf("#4 MIS combined estimator after %d samples is : %f, error: %f (Balance Heuristic) \n",i,MISestimator2,error_MISestimator2);
			printf("#5 MIS combined estimator after %d samples is : %f, error: %f (Power Heuristic, beta=%d)\n",i,MISestimator3,error_MISestimator3,beta);
			printf("#6 MIS combined estimator after %d samples is : %f, error: %f (Power Heuristic, beta=%d)\n",i,MISestimator4,error_MISestimator4,beta_2);

			printf("-------------------\n\n\n");
			report++;
		}
	}
	std::cout << "Press enter to continue ..."; 
	std::cin.get();
	return 0;
}
