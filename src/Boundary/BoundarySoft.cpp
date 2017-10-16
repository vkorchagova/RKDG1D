#include "BoundarySoft.h"

//Конструктор периодических ГУ
BoundarySoft::BoundarySoft(const BaseParams& prm, const Problem& prb) : Boundary(prm, prb)
{
}

//Деструктор периодических ГУ
BoundarySoft::~BoundarySoft()
{
}

void BoundarySoft::FillVirtualCells(const vector<vector<vector<double>>>& SOL) const
{
	int nx = ptrprm->nx;

	SET_NO_CONST_REF(LSOL, SOL[nx]);
	SET_NO_CONST_REF(RSOL, SOL[nx + 1]);
	
	const vector<double>& S0L = SOL[0][0];
	const vector<double>& S0R = SOL[nx - 1][0];

	if (dynamic_cast<const ProblemGas1D*>(ptrprb))
	{
		LSOL = { { S0L[var::r], S0L[var::rvx], S0L[var::rvy], S0L[var::rvz], S0L[var::e] } };
		RSOL = { { S0R[var::r], S0R[var::rvx], S0R[var::rvy], S0R[var::rvz], S0R[var::e] } };

		if (ptrprb->nshape > 1)
		{
			const vector<double>& S1L = SOL[0][1];
			const vector<double>& S1R = SOL[nx - 1][1];

			LSOL.push_back({ 0.0, 0.0, 0.0, 0.0, 0.0 });
			RSOL.push_back({ 0.0, 0.0, 0.0, 0.0, 0.0 });
		}

		if (ptrprb->nshape > 2)
		{
			const vector<double>& S2L = SOL[0][2];
			const vector<double>& S2R = SOL[nx - 1][2];

			LSOL.push_back({ 0.0, 0.0, 0.0, 0.0, 0.0 });
			RSOL.push_back({ 0.0, 0.0, 0.0, 0.0, 0.0 });
		}

		//		cout << "LW, RW !!!" << endl;
	}

	//CopySol(SOL[nx - 1], RSOL);
	//CopySol(SOL[0], LSOL);
}
