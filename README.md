Within-host stochastic viral dynamics code
==========================================

Description
-----------
This repo contains a couple of programmes which may be used to
simulate the early stochastic dynamics of within-host viral infection.
The parameters are set by default to approximate those applicable to
HIV-1 infection dynamics.

File descriptions
-----------------
The most important files are:

* hiv_gillespie.cc: Gillespie SSA simulation code (C++)
* hiv_tauleap_hybrid.cc: Hybrid tau-leaping simulation code (C++)
* deterministic.py: Integrator for deterministic model (Python)
* poissonian.cc: Methods to sample from Poissonian distributions (C++)

Requirements
------------
Both the SSA and tau-leaping programmes use MPI to coordinate the
parallel calculation of trajectories. (Multithreading is adequate for
a single machine, but MPI allows one to easily use a cluster as well.)
So, you'll need some MPI implementation installed, even if you don't
plan to use the parallel feature.

Directions
----------
Edit one of the above files to ensure the parameters are close to what
you want, then use make to generate the executables. Running each
executable without any arguments generates instructions on how to use
it.

