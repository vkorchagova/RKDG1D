#include "TimestepRK2.h"

TimestepRK2::TimestepRK2(const BaseParams& prm, int dimension, int nshape, const Boundary& bnd)
: Timestep(prm, dimension, nshape, bnd)
{
}

TimestepRK2::~TimestepRK2()
{
}

//Шаг расчета методом Рунге-Кутты 2-го порядка
void TimestepRK2::runstep(const vector<vector<vector<double>>>& SOL, \
	vector<vector<vector<double>>>& SOLnew, \
	const double tau, Flux& method, Limiter& lim)
{
	//substep 1
	ApplyBoundary(SOL);
	method.step(SOL, DSOL, 0.5*tau);

	SOLnew = SOL;
	SOLnew += DSOL;

	ApplyBoundary(SOLnew);
	lim.Bound(SOLnew);

	//substep 2
	ApplyBoundary(SOLnew);
	method.step(SOLnew, DSOL, tau);

	SOLnew = SOL; 
	SOLnew += DSOL;

	ApplyBoundary(SOLnew);
	lim.Bound(SOLnew);
}
