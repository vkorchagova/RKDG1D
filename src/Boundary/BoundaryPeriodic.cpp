#include "BoundaryPeriodic.h"

//Конструктор периодических ГУ
BoundaryPeriodic::BoundaryPeriodic(const BaseParams& prm, const Problem& prb) : Boundary(prm, prb)
{
}

//Деструктор периодических ГУ
BoundaryPeriodic::~BoundaryPeriodic()
{
}

void BoundaryPeriodic::FillVirtualCells(const vector<vector<vector<double>>>& SOL) const
{
	int nx = ptrprm->nx;

	SET_NO_CONST_REF(LSOL, SOL[nx]);	
	SET_NO_CONST_REF(RSOL, SOL[nx + 1]);
	
	CopySol(SOL[nx - 1], LSOL);
	CopySol(SOL[0], RSOL);
}
