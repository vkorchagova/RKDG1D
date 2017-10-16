#include "LimiterFinDiff.h"

//КОНСТРУКТОР Лимитера FinDiff
LimiterFinDiff::LimiterFinDiff(const BaseParams& prm, const Problem& prb,\
    const Indicator& ind) : Limiter(prm, prb, ind)
{};

//ДЕСТРУКТОР Лимитера
LimiterFinDiff::~LimiterFinDiff()
{};


void LimiterFinDiff::CalculateBound(const vector<vector<vector<double>>>& SOL, const int cell)
{
	int dim = ptrprb->dim;
    for (int val = 0; val < dim; ++val)
		for (int shape = 1; shape < ptrprb->nshape; ++shape)
			SOLcorr[cell][shape][val] = 0.0;    
}
