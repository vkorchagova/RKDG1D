#include "Timestep.h"

//Конструктор
Timestep::Timestep(const BaseParams& prm, int dimension, int nshape, const Boundary& bnd)
{
	ptrprm = &prm;
	dim = dimension;
	ptrbnd = &bnd;

	
	int nx = ptrprm->nx;

	DSOL.resize(nx + 2);

	for (int cell = 0; cell < (nx + 2); ++cell)
	{
		DSOL[cell].resize(nshape);
		
		for (int shape = 0; shape < nshape; ++shape)
			DSOL[cell][shape].resize(dim);
	}	
}

//Деструктор
Timestep::~Timestep()
{
}

void Timestep::ApplyBoundary(const vector<vector<vector<double>>>& SOL)
{	
	ptrbnd->FillVirtualCells( SOL );
}
