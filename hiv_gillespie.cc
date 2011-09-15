// hiv_gillespie.cc: Implements Gillespie's method to directly integrate a birth/death
// master equation.  This version supports parallel execution of trajectories via MPI.
//
//	X --> 0 (rate d)
//	Y --> 0 (rate a)
//	V --> 0 (rate u)
//
//	0 --> X (rate lambda)
//	
//	X + V --> Y (rate beta)
//
//	Y --> V + Y (rate k)
//	
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cmath>
#include <cstdlib>

#include <mpi.h>

#include <ctime>

#include "poissonian.h"

class Moment {

	public:
		std::string name;
		std::vector <int> exp;
		std::vector <double> mean;
		std::vector <double> sem;

		// Constructor:
		Moment (unsigned int NS, unsigned int Nsamples)
		{
			exp.resize(NS);
			mean.resize(Nsamples, 0.0);
			sem.resize(Nsamples, 0.0);
		}

		// Sampling
		void sample(std::vector<double> X, unsigned int samp)
		{

			double thismoment = 1.0;
			for (unsigned int s=0; s<X.size(); s++) {
				for (unsigned int i=0; i<exp[s]; i++)
					thismoment *= X[s];

			}
			mean[samp] += thismoment;
			sem[samp] += thismoment*thismoment;
		}

		// Normalisation:
		void normalise(unsigned int Npaths)
		{
			for (int s=0; s<mean.size(); s++) {
				mean[s] /= (double)Npaths;
				sem[s] /= (double)Npaths;
				sem[s] = sqrt((sem[s] - mean[s]*mean[s])/(double)Npaths);
			}

		}
};

void sample(std::vector<double> X, std::vector<Moment> & moments, unsigned int samp)
{
	for (unsigned int m=0; m<moments.size(); m++)
		moments[m].sample(X, samp);
}

int main(int argc, char **argv)
{
	using namespace std;

	if (argc < 3) {
		cout << "Usage: " << argv[0] << " <seed> <output file>" << endl;
		exit(0);
	}

	unsigned short cl_seed = (unsigned short)strtol(argv[1], NULL, 10);
	char *ofilename = argv[2];

	// Initialise simulation parameters:

	double T = 30; // Length of simulation

	unsigned int Nsamples = 1001; // Number of samples
	unsigned int Npaths = 400; // Number of independent trajectories
	double samp_dt = T/((double)(Nsamples - 1));


	// Define number and type of reaction processes:

	unsigned int NP = 6; // Number of reaction processes
	unsigned int NS = 3; // Number of species
	unsigned int *R = new unsigned int[2*NS*NP];
	double *c = new double[NP];
	std::vector<double> X(NS, 0.0);

	// Uninfected T cell production
	R[0] = 0; R[1] = 0; R[2] = 0;
	// -->
	c[0] = 1e5; 
	R[3] = 1; R[4] = 0; R[5] = 0;

	// Uninfected T cell death
	R[6] = 1; R[7] = 0; R[8] = 0;
	// -->
	c[1] = 0.1; 
	R[9] = 0; R[10] = 0; R[11] = 0;

	// Infected T cell death
	R[12] = 0; R[13] = 1; R[14] = 0;
	// -->
	c[2] = 0.5; 
	R[15] = 0; R[16] = 0; R[17] = 0;

	// Virion death
	R[18] = 0; R[19] = 0; R[20] = 1;
	// -->
	c[3] = 5.0; 
	R[21] = 0; R[22] = 0; R[23] = 0;

	// T cell infection
	R[24] = 1; R[25] = 0; R[26] = 1;
	// -->
	c[4] = 2e-7; 
	R[27] = 0; R[28] = 1; R[29] = 0;

	// Virion production
	R[30] = 0; R[31] = 1; R[32] = 0;
	// -->
	c[5] = 100.0; 
	R[33] = 0; R[34] = 1; R[35] = 1;

	// Set up output moments:

	std::vector<Moment> moments;
	Moment thismoment(NS, Nsamples);

	thismoment.name = "x";
	thismoment.exp[0] = 1; thismoment.exp[1] = 0; thismoment.exp[2] = 0;
	moments.push_back(thismoment);

	thismoment.name = "y";
	thismoment.exp[0] = 0; thismoment.exp[1] = 1; thismoment.exp[2] = 0;
	moments.push_back(thismoment);

	thismoment.name = "v";
	thismoment.exp[0] = 0; thismoment.exp[1] = 0; thismoment.exp[2] = 1;
	moments.push_back(thismoment);

	thismoment.name = "x2";
	thismoment.exp[0] = 2; thismoment.exp[1] = 0; thismoment.exp[2] = 0;
	moments.push_back(thismoment);

	thismoment.name = "y2";
	thismoment.exp[0] = 0; thismoment.exp[1] = 2; thismoment.exp[2] = 0;
	moments.push_back(thismoment);

	thismoment.name = "v2";
	thismoment.exp[0] = 0; thismoment.exp[1] = 0; thismoment.exp[2] = 2;
	moments.push_back(thismoment);

	thismoment.name = "xy";
	thismoment.exp[0] = 1; thismoment.exp[1] = 1; thismoment.exp[2] = 0;
	moments.push_back(thismoment);

	thismoment.name = "xv";
	thismoment.exp[0] = 1; thismoment.exp[1] = 0; thismoment.exp[2] = 1;
	moments.push_back(thismoment);

	thismoment.name = "yv";
	thismoment.exp[0] = 0; thismoment.exp[1] = 1; thismoment.exp[2] = 1;
	moments.push_back(thismoment);

	// Keep track of total number of reactions performed:
	double Nreactions = 0;


	// Declare temporary variables necessary for algorithm:

	double *a = new double[NP]; 
	double a0, asum, h;

	// MPI Initialisation boilerplate:

	MPI::Init(argc, argv);
	int mpi_rank = MPI::COMM_WORLD.Get_rank();
	int mpi_size = MPI::COMM_WORLD.Get_size();


	// Determine paths for this MPI process to integrate:
	unsigned int chunk_size = Npaths / mpi_size;
	if (Npaths % mpi_size > 0) {
		chunk_size++;
		Npaths = chunk_size*mpi_size;
	}

	// Initialise PRNG:

	unsigned short buf[3];
	buf[0] = 4253;
	buf[1] = cl_seed;
	buf[2] = mpi_rank;


	// Start timer:

	time_t start_time, end_time;
	time(&start_time);


	// Loop over paths in ensemble:
	for (unsigned int path=0; path<chunk_size; path++) {

		// Report path number:
		cout << "Rank " << mpi_rank << ": Starting path " << path+1 << " of " << chunk_size << endl;

		// Initialise sampling index:
		unsigned int samp = 0;

		// Set initial time:
		double t = 0;

		// Set Poissonian initial conditions for subpopulations:
		X[0] = poissonian(1e6, buf);
		X[1] = poissonian(0.0, buf);
		X[2] = poissonian(1e2, buf);


		// Sample initial state:
		sample(X, moments, samp++);


		////////////////////////////////////////
		// SIMULATION LOOP

		while (true) {

			// Step 1: calculate a[]

			a0 = 0.0;
			for (unsigned int p=0; p<NP; p++) {
				h = 1.0;
				for (unsigned int s=0; s<NS; s++) {
					for (unsigned int j=0; j<R[s + p*2*NS]; j++)
						h *= (X[s]-(double)j);
				}
				a[p] = h*c[p];
				a0 += a[p];
			}

			// Step 2: increment t and find mu

			double r1 = erand48(buf);
			double r2 = erand48(buf);

			t += -(1.0/a0)*log(r1);
			asum = 0.0;
			unsigned int mu = 0;
			for (unsigned int p=0; p<NP; p++) {
				asum += a[p];
				if (asum >= r2*a0) {
					mu = p;
					break;
				}
			}

			// Step 3: Sample result at regular intervals

			while ((samp_dt*(double)samp < t) && (samp_dt*(double)samp <= T))
				sample(X, moments, samp++);

			// Step 4: Break if necessary:
			if (t>T)
				break;


			// Step 5: Implement chosen reaction

			for (unsigned int s=0; s<NS; s++) {
				X[s] -= (double)R[s + mu*2*NS];
				X[s] += (double)R[s + NS + mu*2*NS];
			}

			// Step 6: Increment reaction counter:

			Nreactions++;

		}

		////////////////////////////////////////

	}

	if (mpi_rank > 0) {

		cout << "Rank " << mpi_rank << ": Sending results" << endl;
		for (int m=0; m<moments.size(); m++) {
			MPI::COMM_WORLD.Send(&(moments[m].mean[0]), Nsamples, MPI::DOUBLE, 0, 2*m);
			MPI::COMM_WORLD.Send(&(moments[m].sem[0]), Nsamples, MPI::DOUBLE, 0, 2*m+1);
		}
		MPI::COMM_WORLD.Send(&Nreactions, 1, MPI::DOUBLE, 0, 2*moments.size());

	} else {
		
		cout << "Rank " << mpi_rank << ": Collating results" << endl;
		double *mean_recv = new double[Nsamples];
		double *sem_recv = new double[Nsamples];
		double Nreactions_recv;

		for (unsigned int recv_rank=1; recv_rank<mpi_size; recv_rank++) {

			for (int m=0; m<moments.size(); m++) {
				MPI::COMM_WORLD.Recv(mean_recv, Nsamples, MPI::DOUBLE, recv_rank, 2*m);
				MPI::COMM_WORLD.Recv(sem_recv, Nsamples, MPI::DOUBLE, recv_rank, 2*m+1);

				for (unsigned int s=0; s<Nsamples; s++) {
					moments[m].mean[s] += mean_recv[s];
					moments[m].sem[s] += sem_recv[s];
				}
			}
			MPI::COMM_WORLD.Recv(&Nreactions_recv, 1, MPI::DOUBLE, recv_rank, 2*moments.size());
			Nreactions += Nreactions_recv;

			cout << "Rank " << mpi_rank << ": Recieved results from " << recv_rank << endl;
		}

		delete [] mean_recv;
	   	delete [] sem_recv;


		// Stop timer:
		time(&end_time);

		// Normalise output moments:
		for (unsigned int m=0; m<moments.size(); m++)
			moments[m].normalise(Npaths);
		Nreactions /= (double)Npaths;


		// Output to file:
		ofstream ofile;
		ofile.open(ofilename);

		ofile << "# output from hiv_gillespie" << endl
			<< "#" << endl
			<< "# Npaths = " << Npaths << endl
			<< "# mpi_size = " << mpi_size << endl
			<< "#" << endl
			<< "# Simulation took "
			<< (unsigned int)(difftime(end_time, start_time))/60
			<< " minutes and "
			<< (unsigned int)(difftime(end_time, start_time)) % 60
			<< " seconds to run." << endl
			<< "#" << endl
			<< "# Mean number of reactions simulated per path: " << Nreactions << endl
			<< "#" << endl
			<< endl;

		ofile << "t ";
		for (int m=0; m<moments.size(); m++) {
			ofile << moments[m].name << " ";
			ofile << moments[m].name << "_sem ";
		}
		ofile << endl;

		for (int s=0; s< Nsamples; s++) {
			ofile << ((double)s)*samp_dt << ' ';
			for (int m=0; m<moments.size(); m++) {
				ofile << moments[m].mean[s] << " ";
				ofile << moments[m].sem[s] << " ";
			}
			ofile << endl;
		}

		ofile.close();

	}

	// Free dynamically allocated memory:
	delete [] R;
	delete [] c;
	delete [] a;

	// MPI cleanup boilerplate:
	MPI::Finalize();

	return 0;
}
