#include "Problem.h"


Problem::Problem(const BaseParams& prm, int dimension, int nshapefunctions, vector<std::function<double(const double)>> initv)
{
	ptrprm = &prm;
	dim = dimension;
	nshape = nshapefunctions;
	init = initv;

	int nx = ptrprm->nx;

	//Выделяем память под все необходимые переменные
	FLUX.resize(nx + 2); LFLUX.resize(nx + 2); RFLUX.resize(nx + 2); // +2 - для фиктивных ячеек

	for (int cell = 0; cell < (nx + 2); ++cell)
	{
		FLUX[cell].resize(dim); LFLUX[cell].resize(dim); RFLUX[cell].resize(dim);
	}

	shapefunc[0] = [](double pt){ return 1.0; };
	shapefunc[1] = [](double pt){ return pt; };
	shapefunc[2] = [](double pt){ return (1.5 * sqr(pt) - 0.5); };

	shapefuncnorm2[0] = 1.0;
	shapefuncnorm2[1] = 1.0/3.0;
	shapefuncnorm2[2] = 1.0/5.0;

}


Problem::~Problem()
{
}






//Запонение вектора потоков для всех ячеек
/*
void Problem::convFlux(const vector<vector<double>>& UU)
{
	for (int cell = 0; cell < (ptrprm->nx + 2); ++cell)
	{
		getFlux(UU[cell], FLUX[cell]);
	}
}
*/

void Problem::convFlux(const vector<vector<vector<double>>>& SOL)
{
	vector<double> UlimL(dim), UlimR(dim), Ucenter(dim);
	
	for (int cell = 0; cell < (ptrprm->nx + 2); ++cell)
	{
		const vector<vector<double>>& sol = SOL[cell];

		for (var val = (var)0; val < dim; val = (var)(val + 1))
		{
			UlimL[val]   = side_val(sol, val, side::left);
			UlimR[val]   = side_val(sol, val, side::right);
			Ucenter[val] = sol[0][val];
			
			if (nshape > 2)
				Ucenter[val] -= 0.5*sol[2][val];

		}

		getFlux(Ucenter, FLUX[cell]);
		getFlux(UlimL,   LFLUX[cell]);
		getFlux(UlimR,   RFLUX[cell]);
	}
}
