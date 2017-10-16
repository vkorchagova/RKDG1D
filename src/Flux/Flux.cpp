#include "Flux.h"


Flux::Flux(const BaseParams& prm, Problem& prb) 
{
	ptrprm = &prm; 
	ptrprb = &prb;

	ptrprb_toGas = dynamic_cast<ProblemGas1D*>(ptrprb);
	ptrprb_toMHD = dynamic_cast<ProblemMHD1D*>(ptrprb);
}


Flux::~Flux()
{
}



void Flux::setlocsolflux(const vector<vector<vector<double>>>& SOL, const int cell)
{
	int nx = ptrprm->nx;

	ptr_mysol = &(SOL[cell]);
	ptr_leftsol = (cell > 0) ? &(SOL[cell - 1]) : &(SOL[nx]);
	ptr_rightsol = (cell < nx - 1) ? &(SOL[cell + 1]) : &(SOL[nx + 1]);


	ptr_leftfluxR = (cell > 0) ? &(ptrprb->RFLUX[cell - 1]) : &(ptrprb->RFLUX[nx]);
	ptr_myflux = &(ptrprb->FLUX[cell]);
	ptr_myfluxL = &(ptrprb->LFLUX[cell]);
	ptr_myfluxR = &(ptrprb->RFLUX[cell]);
	ptr_rightfluxL = (cell < nx - 1) ? &(ptrprb->LFLUX[cell + 1]) : &(ptrprb->LFLUX[nx + 1]);
}




