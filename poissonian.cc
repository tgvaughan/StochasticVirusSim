// poisson.cc - Implements routines for generating Poisson-distributed random
// variables.  Note: values are returned as doubles rather than unsigned ints
// which would perhaps be more appropriate but I think is just finicky.

#include <stdlib.h>
#include <cmath>

#include "poissonian.h"

#define PI 3.141592653589793238


// If you believe NR, returns the value ln[Gamma(xx)] for xx>0. Pure black magick.
double gammln(double xx)
{
	int j;
	double x,y,tmp,ser;
	static const double cof[6]={76.18009172947146, -86.50532032941677,
		24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
		-0.5395239384953e-5};

	y = x = xx;

	tmp = x + 5.5;
	tmp -= (x + 0.5)*log(tmp);

	ser=1.000000000190015;

	for (j=0; j<6; j++)
		ser += cof[j]/++y;

	return -tmp + log(2.5066282746310005*ser/x);
}

// Rejection method from NR, apparently good for lambda>=12
double poissonian_reject(double lambda, unsigned short *buf)
{
	double sq = sqrt(2.0*lambda);
	double alxm = log(lambda);
	double g = lambda*alxm - gammln(lambda+1.0);
	double em, t, y;

	do {
		do {
			y = tan(PI*erand48(buf));
			em = sq*y + lambda;
		} while (em < 0.0);

		em = floor(em);
		t = 0.9*(1.0 + y*y)*exp(em*alxm - gammln(em + 1.0)-g);

	} while (erand48(buf) > t);

	return em;
}

// Direct method due to Knuth. Only efficient for small lambda. Bad Knuth.
double poissonian_knuth(double lambda, unsigned short *buf)
{
	double L = exp(-lambda);
	double p;
	unsigned int k;
	
	for (k=0, p=1; p >= L; k++)
		p = p*erand48(buf);

	return k-1;
}

// Meta-function selects direct or rejection method according to lambda.
double poissonian(double lambda, unsigned short *buf)
{
	if (lambda < 12.0)
		return poissonian_knuth(lambda, buf);

	return poissonian_reject(lambda, buf);
}
