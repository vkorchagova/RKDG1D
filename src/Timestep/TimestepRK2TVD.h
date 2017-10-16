#ifndef TIMESTEPRK2TVD_H_
#define TIMESTEPRK2TVD_H_

#include "Timestep.h"

#define RK2TVD TimestepRK2TVD

class TimestepRK2TVD :
	public Timestep
{
protected:
	//¬екторы дл€ хранени€ промежуточных значений:
	vector<vector<vector<double>>> DSOLtemp;

public:
	TimestepRK2TVD(const BaseParams& prm, int dimension, int nshape, const Boundary& bnd);
	~TimestepRK2TVD();

	void runstep(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& SOLnew, \
		const double tau, Flux& method, Limiter& lim);
};

#endif
