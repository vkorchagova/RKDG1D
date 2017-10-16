#include "TimestepRK2TVD.h"

TimestepRK2TVD::TimestepRK2TVD(const BaseParams& prm, int dimension, int nshape, const Boundary& bnd)
	: Timestep(prm, dimension, nshape, bnd)
{
	int nx = ptrprm->nx;

	DSOLtemp.resize(nx + 2);

	for (int cell = 0; cell < (nx + 2); ++cell)
	{
		DSOLtemp[cell].resize(nshape);

		for (int shape = 0; shape < nshape; ++shape)
			DSOLtemp[cell][shape].resize(dim);
	}
}

TimestepRK2TVD::~TimestepRK2TVD()
{
}

//Шаг расчета методом Рунге-Кутты 2-го порядка (схема TVD)
void TimestepRK2TVD::runstep(const vector<vector<vector<double>>>& SOL, \
	vector<vector<vector<double>>>& SOLnew, \
	const double tau, Flux& method, Limiter& lim)
{
	//substep 1
	ApplyBoundary(SOL);
	method.step(SOL, DSOL, tau);

	SOLnew = SOL;
	SOLnew += DSOL;
		
	DSOLtemp = DSOL;

	ApplyBoundary(SOLnew);
	lim.Bound(SOLnew);

	//substep 2
	ApplyBoundary(SOLnew);
	method.step(SOLnew, DSOL, tau);

	SOLnew = SOL;
	DSOL += DSOLtemp;
	DSOL *= 0.5;
	SOLnew += DSOL;

	ApplyBoundary(SOLnew);
	lim.Bound(SOLnew);
}