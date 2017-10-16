#include "TimestepRK3TVD.h"

TimestepRK3TVD::TimestepRK3TVD(const BaseParams& prm, int dimension, int nshape, const Boundary& bnd)
: Timestep(prm, dimension, nshape, bnd)
{
	int nx = ptrprm->nx;


	SOLtemp1.resize(nx + 2);
	SOLtemp2.resize(nx + 2);

	for (int cell = 0; cell < (nx + 2); ++cell)
	{
		SOLtemp1[cell].resize(nshape);
		SOLtemp2[cell].resize(nshape);

		for (int shape = 0; shape < nshape; ++shape)
		{
			SOLtemp1[cell][shape].resize(dim);
			SOLtemp2[cell][shape].resize(dim);
		}
	}
}

TimestepRK3TVD::~TimestepRK3TVD()
{
}

//Шаг расчета методом Рунге-Кутты 3-го порядка (схема TVD)
void TimestepRK3TVD::runstep(const vector<vector<vector<double>>>& SOL, \
	vector<vector<vector<double>>>& SOLnew, \
	const double tau, Flux& method, Limiter& lim)
{
	//substep 1
	ApplyBoundary(SOL);
	method.step(SOL, DSOL, tau);

	SOLtemp1 = SOL;
	SOLtemp1 += DSOL;
	
	ApplyBoundary(SOLtemp1);
	lim.Bound(SOLtemp1);

	//substep 2
	ApplyBoundary(SOLtemp1);
	method.step(SOLtemp1, DSOL, 0.25*tau);

	SOLtemp2 = SOL*0.75;
	SOLtemp2 += SOLtemp1*0.25;
	SOLtemp2 += DSOL;	

	ApplyBoundary(SOLtemp2);
	lim.Bound(SOLtemp2);


	//substep 3
	ApplyBoundary(SOLtemp2);
	method.step(SOLtemp2, DSOL, 0.66666666666666666666666*tau);

	SOLnew = SOL*0.3333333333333333333;
	SOLnew += SOLtemp2*0.66666666666666666666;
	SOLnew += DSOL;	

	ApplyBoundary(SOLnew);
	lim.Bound(SOLnew);
}