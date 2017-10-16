#include "LimiterWENO.h"
#include "Integrator.h"


//КОНСТРУКТОР Лимитера WENO
LimiterWENO::LimiterWENO(const BaseParams& prm, const Problem& prb,\
    const Indicator& ind, double degree) : Limiter(prm, prb, ind)
{
    wg = degree;
}

//ДЕСТРУКТОР Лимитера
LimiterWENO::~LimiterWENO()
{}


void LimiterWENO::CalculateBound(const vector<vector<vector<double>>>& SOL, const int cell)
{
	int nx = ptrprm->nx;
	double h = ptrprm->h;
	int dim = ptrprb->dim;

	const vector<double>& leftU = (cell > 0) ? SOL[cell - 1][0] : SOL[nx][0];
	const vector<double>& rightU = (cell < nx - 1) ? SOL[cell + 1][0] : SOL[nx + 1][0];

	const vector<double>& leftleftU = (cell > 1) ? SOL[cell - 2][0] : SOL[nx][0];
	const vector<double>& rightrightU = (cell < nx - 2) ? SOL[cell + 2][0] : SOL[nx + 1][0];

	const vector<vector<double>> linwei = \
	{{ 0.5625, 0.375, 0.0625 },
	{ 0.374303, 0.475, 0.150697 },
	{ 0.150697, 0.475, 0.374303 },
	{ 0.0625, 0.375, 0.5625 } };


	auto& shapes = ptrprb->shapefunc;
	auto& shapesnrms2 = ptrprb->shapefuncnorm2;

	for (int val = 0; val < dim; ++val)
	{
		const double& yll = leftleftU[val];
		const double& yl  = leftU[val];
		const double& y   = SOL[cell][0][val];
		const double& yr  = rightU[val];
		const double& yrr = rightrightU[val];


		double beta[3];
		double dif[3] = { backbackdiff(yll, yl, y, h), centrdiff(yl, yr, h), forwforwdiff(y, yr, yrr, h) };
		double difdif[3] = { 0.5*diff2(yll, yl, y, h), 0.5*diff2(yl, y, yr, h), 0.5*diff2(y, yr, yrr, h) };
		for (int wei = 0; wei < 3; wei++)
			beta[wei] = sqr(h*dif[wei]) + 13.0 *sqr(h*h*difdif[wei]) / 3.0;

		for (int shape = 1; shape < ptrprb->nshape; shape++)
			SOLcorr[cell][shape][val] = 0.0;

		for (int gp = 0; gp < 4; gp++)
		{
			double xi = 0.5 * h * pos[Lobatto][3][gp];
			double p[3];
			double nonlinwei[3];

			double sumwei = 0.0;

			for (int wei = 0; wei < 3; wei++)
			{
				p[wei] = y + xi*dif[wei] + sqr(xi)*difdif[wei];
				nonlinwei[wei] = linwei[gp][wei] / pow(weps + beta[wei], wg);
				sumwei += nonlinwei[wei];
			}

			double combVal = 0.0;
			for (int wei = 0; wei < 3; wei++)
				combVal += nonlinwei[wei] * p[wei] / sumwei;

			for (int shape = 1; shape < ptrprb->nshape; shape++)
				SOLcorr[cell][shape][val] += wgt[Lobatto][3][gp] *  0.5*combVal*shapes[shape](pos[Lobatto][3][gp]) / shapesnrms2[shape];
		}
	}
}
