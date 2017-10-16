#ifndef TIMESTEPRK2_H_
#define TIMESTEPRK2_H_

#include "Timestep.h"

#define RK2 TimestepRK2

class TimestepRK2 :
	public Timestep
{
public:
	TimestepRK2(const BaseParams& prm, int dimension, int nshape, const Boundary& bnd);
	~TimestepRK2();

	void runstep(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& SOLnew, \
		const double tau, Flux& method, Limiter& lim);
};

#endif

