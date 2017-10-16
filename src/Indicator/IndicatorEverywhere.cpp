#include "IndicatorEverywhere.h"

//КОНСТРУКТОР индикатора "везде и всюду"
IndicatorEverywhere::IndicatorEverywhere(const BaseParams& prm, const Problem& prb) : Indicator(prm, prb)
{
}//IndicatorEverywhere::IndicatorEverywhere

//ДЕСТРУКТОР индикатора "везде и всюду"
IndicatorEverywhere::~IndicatorEverywhere()
{};

//Вычисление индикаторной функции
void IndicatorEverywhere::calc_indicator(const vector<vector<vector<double>>>& SOL, \
	vector<double>& Ind) const
{
	int nx = ptrprm->nx;
	double h = ptrprm->h;

	for (int cell = 0; cell < nx; ++cell)
		Ind[cell] = 2.0;
}
