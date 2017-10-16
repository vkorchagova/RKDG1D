#include "ProblemTransfer1D.h"

const int dimTransfer1D = 1;

ProblemTransfer1D::ProblemTransfer1D(const BaseParams& prm, int order, \
	vector<std::function<double(const double)>> initv, \
	std::function<double(const double)> vel) : Problem(prm, dimTransfer1D, order, initv)
{
	physflux = vel;
}


ProblemTransfer1D::~ProblemTransfer1D()
{
}


void ProblemTransfer1D::getFlux(const vector<double>& U, vector<double>& Flux) const
{
	Flux[var::U] = physflux(U[var::U]);
}


//Значения всех переменных на левой и правой границах ячеек
inline double ProblemTransfer1D::side_val(const vector<vector<double>>& sol, var q, side sd) const
{
	double sgn = (sd == side::left) ? -1.0 : 1.0;
	double res = 0.0;

	switch (q)
	{
	case var::U:
		res = sol[0][q] + sgn * sol[1][q];
		if (nshape>2)
			res += sol[2][q];
		break;
	}//switch
	return res;
}



//Начальные условия
double ProblemTransfer1D::initial_var(const double x, const var q, const int cft) const
{
	const int add = cft*dim; //(cft == 0) ? 0 : dim;
	double res = 0.0;

	//Для удобства ссылки:
	// = ptrprm->basic.init;

	switch (q)
	{
		//Условия для U задаем непосредственно
	case var::U:
		res = (init[0 + add])(x);
		break;

	}//switch q

	return res;
}//initial_var

