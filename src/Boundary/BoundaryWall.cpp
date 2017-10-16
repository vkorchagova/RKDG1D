#include "BoundaryWall.h"

//Конструктор ГУ твердой стенки
BoundaryWall::BoundaryWall(const BaseParams& prm, const Problem& prb) : Boundary(prm, prb)
{
}

//Деструктор ГУ твердой стенки
BoundaryWall::~BoundaryWall()
{
}

void BoundaryWall::FillVirtualCells(const vector<vector<vector<double>>>& SOL) const
{
	int nx = ptrprm->nx;

	SET_NO_CONST_REF(LSOL, SOL[nx]);
	SET_NO_CONST_REF(RSOL, SOL[nx + 1]);
	
	const vector<double>& S0L = SOL[0][0];
	const vector<double>& S0R = SOL[nx - 1][0];

	if (dynamic_cast<const ProblemGas1D*>(ptrprb))
	{
		LSOL = { { S0L[var::r], -S0L[var::rvx], S0L[var::rvy], S0L[var::rvz], S0L[var::e] } };
		RSOL = { { S0R[var::r], -S0R[var::rvx], S0R[var::rvy], S0R[var::rvz], S0R[var::e] }	};

		if (ptrprb->nshape > 1)
		{
			const vector<double>& S1L = SOL[0][1];
			const vector<double>& S1R = SOL[nx - 1][1];
			
			LSOL.push_back( { -S1L[var::r], S1L[var::rvx], -S1L[var::rvy], -S1L[var::rvz], -S1L[var::e] } );
			RSOL.push_back(	{ -S1R[var::r], S1R[var::rvx], -S1R[var::rvy], -S1R[var::rvz], -S1R[var::e] } );
		}
		
		if (ptrprb->nshape > 2)
		{
			const vector<double>& S2L = SOL[0][2];
			const vector<double>& S2R = SOL[nx - 1][2];

			LSOL.push_back(	{ 0.0, 0.0, 0.0, 0.0, 0.0 } );
			RSOL.push_back( { 0.0, 0.0, 0.0, 0.0, 0.0 } );
		}

//		cout << "LW, RW !!!" << endl;
	}

	if (dynamic_cast<const ProblemMHD1D*>(ptrprb))
	{
		LSOL = { { S0L[var::r], -S0L[var::rvx], S0L[var::rvy], S0L[var::rvz], S0L[var::e], S0L[var::Hy], S0L[var::Hz]} };
		RSOL = { { S0R[var::r], -S0R[var::rvx], S0R[var::rvy], S0R[var::rvz], S0R[var::e], S0R[var::Hy], S0R[var::Hz]} };

		if (ptrprb->nshape > 1)
		{
			const vector<double>& S1L = SOL[0][1];
			const vector<double>& S1R = SOL[nx - 1][1];

			LSOL.push_back({ -S1L[var::r], S1L[var::rvx], -S1L[var::rvy], -S1L[var::rvz], -S1L[var::e], -S1L[var::Hy], -S1L[var::Hz] });
			RSOL.push_back({ -S1R[var::r], S1R[var::rvx], -S1R[var::rvy], -S1R[var::rvz], -S1R[var::e], -S1R[var::Hy], -S1R[var::Hz] });
		}

		if (ptrprb->nshape > 2)
		{
			const vector<double>& S2L = SOL[0][2];
			const vector<double>& S2R = SOL[nx - 1][2];

			LSOL.push_back({ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });
			RSOL.push_back({ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });
		}

		//		cout << "LW, RW !!!" << endl;
	}

	if (dynamic_cast<const ProblemTransfer1D*>(ptrprb))
	{
		LSOL = { { SOL[0][0][var::U] } };
		RSOL = { { SOL[nx - 1][0][var::U] } };

		if (ptrprb->nshape > 1)
		{
			LSOL.push_back({ -SOL[0][1][var::U] });
			RSOL.push_back({ -SOL[nx-1][1][var::U] });
		}
		
		if (ptrprb->nshape > 2)
		{
			LSOL.push_back({ 0.0 });
			RSOL.push_back({ 0.0 });
		}
		
		cout << "LW, RW !!!" << endl;
	}
}
