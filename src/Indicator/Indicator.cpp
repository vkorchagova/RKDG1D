#include "Indicator.h"



Indicator::Indicator(const BaseParams& prm, const Problem& prb)
{
	ptrprm = &prm;
	ptrprb = &prb;
}

Indicator::~Indicator()
{
}
