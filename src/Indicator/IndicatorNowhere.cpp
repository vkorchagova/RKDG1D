#include "IndicatorNowhere.h"

//КОНСТРУКТОР индикатора "никогда"
IndicatorNowhere::IndicatorNowhere(const BaseParams& prm, const Problem& prb) : Indicator(prm, prb)
{
}//IndicatorEverywhere::IndicatorNowhere

//ДЕСТРУКТОР индикатора "никогда"
IndicatorNowhere::~IndicatorNowhere()
{};

//Вычисление индикаторной функции
void IndicatorNowhere::calc_indicator(const vector<vector<vector<double>>>& SOL, \
	vector<double>& Ind) const
{
	int nx = ptrprm->nx;
	double h = ptrprm->h;

	for (int cell = 0; cell < nx; ++cell)
		Ind[cell] = 0.0;
}
