#ifndef TIMESTEPEULER_H_
#define TIMESTEPEULER_H_

#include "Timestep.h"

#define Euler TimestepEuler

class TimestepEuler :
	public Timestep
{
public:
	TimestepEuler(const BaseParams& prm, int dimension, int nshape, const Boundary& bnd);
	~TimestepEuler();

	void runstep(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& SOLnew, \
		const double tau, Flux& method, Limiter& lim);
};

#endif