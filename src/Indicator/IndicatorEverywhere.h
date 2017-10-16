#ifndef INDICATOREVERYWHERE_H_
#define INDICATOREVERYWHERE_H_

#include "Indicator.h"
#include "defs.h"

#define Everywhere IndicatorEverywhere 

class IndicatorEverywhere :
	public Indicator
{
public:
	IndicatorEverywhere(const BaseParams& prm, const Problem& prb);
	~IndicatorEverywhere();

	//Вычисление индикаторной функции
	void calc_indicator(const vector<vector<vector<double>>>& SOL, \
		vector<double>& Ind) const;
};

#endif