#ifndef TIMESTEPRK3TVD_H_
#define TIMESTEPRK3TVD_H_

#include "Timestep.h"

#define RK3TVD TimestepRK3TVD

class TimestepRK3TVD :
	public Timestep
{
protected:
	//¬екторы дл€ хранени€ промежуточных значений:
	vector<vector<vector<double>>> SOLtemp1;
	vector<vector<vector<double>>> SOLtemp2;

public:
	TimestepRK3TVD(const BaseParams& prm, int dimension, int nshape, const Boundary& bnd);
	~TimestepRK3TVD();

	void runstep(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& SOLnew, \
		const double tau, Flux& method, Limiter& lim);
};

#endif
