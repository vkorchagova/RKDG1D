#include "IndicatorHarten.h"



//КОНСТРУКТОР индикатора Хартена
IndicatorHarten::IndicatorHarten(const BaseParams& prm, const Problem& prb, const var val, const double alpha) : Indicator(prm, prb)
{
	sens = val;
	mult = alpha;
}//IndicatorHarten::IndicatorHarten

//ДЕСТРУКТОР индикатора Хартена
IndicatorHarten::~IndicatorHarten()
{};

//Вычисление индикаторной функции
void IndicatorHarten::calc_indicator(const vector<vector<double>>& UU, \
	const vector<vector<double>>& VV, \
	vector<double>& Ind) const
{
	bool a1, a2, a3, a4;
	
	int nx = ptrprm->nx;
	double h = ptrprm->h;
	
	//Временная переенная: предварительно рассчитанное значение индикаторной функции
	vector<double> preind(nx);

	for (int cell = 1; cell < nx - 1; ++cell)
	{
		double x1 = (UU[cell - 1][sens] + 2.0*VV[cell - 1][sens] - UU[cell][sens]);
		double x2 = (UU[cell + 1][sens] - 2.0*VV[cell + 1][sens] - UU[cell][sens]);
		
		a1 = x1* \
			 x2 < -1e-7;
		a2 = (fabs(VV[cell][sens]) > mult*fabs(VV[cell - 1][sens])) || \
			 (mult*fabs(VV[cell][sens]) < fabs(VV[cell - 1][sens]));
		a3 = (fabs(VV[cell + 1][sens]) > mult*fabs(VV[cell][sens])) || \
			 (mult*fabs(VV[cell + 1][sens]) < fabs(VV[cell][sens]));
		a4 = a1 && (a2 || a3);
		
		preind[cell] = (a4) ? 1.0 : 0.0;
	}

	preind[0] = 1.0;
	preind[nx - 1] = 1.0;

	for (int cell = 1; cell < nx - 1; ++cell)
		Ind[cell] = 2.0*(preind[cell - 1] + preind[cell] + preind[cell + 1]);

	Ind[0] = 1.0;
	Ind[nx - 1] = 1.0;
}