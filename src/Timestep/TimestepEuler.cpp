#include "TimestepEuler.h"

TimestepEuler::TimestepEuler(const BaseParams& prm, int dimension, int nshape, const Boundary& bnd)
	: Timestep(prm, dimension, nshape, bnd)
{
}

TimestepEuler::~TimestepEuler()
{
}

//Ўаг расчета методом Ёйлера
void TimestepEuler::runstep(const vector<vector<vector<double>>>& SOL, \
	vector<vector<vector<double>>>& SOLnew, \
	const double tau, Flux& method, Limiter& lim)
{
	ApplyBoundary(SOL);
	method.step(SOL, DSOL, tau);

	SOLnew = SOL;
	SOLnew += DSOL;
	
	ApplyBoundary(SOLnew);
	lim.Bound(SOLnew);
}