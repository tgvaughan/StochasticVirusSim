#!/usr/bin/python

# Simple ODE integrator for virus problem.

import sys
from numpy import *

class params:
	# My parameters:
	#lam = 2.5e8
	#beta = 5e-13
	#k = 1e3
	#d = 1e-3
	#a = 1.0
	#u = 3.0

	# Nowak and May's parameters:
	#lam = 1e5
	#beta = 2e-7
	#k = 1e2
	#d = 0.1
	#a = 0.5
	#u = 5.0

	# My other parameters:
	V = 2.5e3
	lam = 1e5*V
	beta = 2e-9/V
	k = 1e3
	d = 0.001
	a = 1.0
	u = 3.0

def ode(a, p):

	x = a[0]
	y = a[1]
	v = a[2]

	dxdt = p.lam - p.beta*x*v - p.d*x
	dydt = p.beta*x*v - p.a*y
	dvdt = p.k*y - p.beta*x*v - p.u*v

	#dxdt = 0.
	#dydt = p.beta*x*v
	#dvdt = p.k*y

	return array((dxdt, dydt, dvdt))


def sistep(a0, f, Niter, dt, params):

	a = a0

	for iter in range(Niter):
		a = a0 + 0.5*ode(a, params)*dt
	
	a = 2.*a - a0

	return a


#### Main

T = 1000
dt = 0.001
Nsteps = int(T/dt + 1)
Nsamples = 1001

Niter = 3

a_init = array((params.lam/params.d,0.,1.))

# Regular sampling for log plot:
samptimes = exp(linspace(log(0.1),log(T),Nsamples))

print "t x y v"

sidx = 0
a = a_init
for tidx in range(Nsteps):
	a = sistep(a, ode, Niter, dt, params)

	while sidx < Nsamples and samptimes[sidx] < tidx*dt:
		print "%g %g %g %g" % (tidx*dt, a[0], a[1], a[2])
		sidx += 1
