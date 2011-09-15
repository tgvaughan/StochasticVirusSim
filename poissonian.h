// poisson.h - Declares the function poisson() which returns a random
// number (cast as a double) belonging to a Poissonian distribution.

double gammln(double xx);
double poissonian_knuth(double lambda, unsigned short *buf);
double poissonian_reject(double lambda, unsigned short *buf);
double poissonian(double lambda, unsigned short *buf);
