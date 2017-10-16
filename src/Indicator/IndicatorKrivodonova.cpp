#include "IndicatorKrivodonova.h"

#include "Problem.h"
#include "ProblemGas1D.h"
#include "ProblemTransfer1D.h"
#include "ProblemMHD1D.h"

//КОНСТРУКТОР индикатора Криводоновой
IndicatorKrivodonova::IndicatorKrivodonova(const BaseParams& prm, \
	const Problem& prb, const var val) : Indicator(prm, prb)
{    
	sens = val;
}//IndicatorKrivodonova::IndicatorKrivodonova

//ДЕСТРУКТОР индикатора Криводоновой
IndicatorKrivodonova::~IndicatorKrivodonova()
{};

//Вычисление индикаторной функции
void IndicatorKrivodonova::calc_indicator(const vector<vector<vector<double>>>& SOL, \
	vector<double>& Ind) const
{
    double zn;

	int nx = ptrprm->nx;
	double h = ptrprm->h;


    for (int cell = 1; cell < nx - 1; ++cell)
    {

		Ind[cell] = 0.0;
		
		int ord = ptrprb->nshape - 1;
		zn = 0.5 * pow(h, 0.5*(ord+1)) * (SOL[cell][0][sens] + 0.000001);

		if (dynamic_cast<const ProblemTransfer1D*>(ptrprb))
			Ind[cell] += ptrprb->side_val(SOL[cell], sens, side::left) - ptrprb->side_val(SOL[cell-1], sens, side::right);

		if (dynamic_cast<const ProblemGas1D*>(ptrprb))
		{

			
			double vleftr = ptrprb->side_val(SOL[cell - 1], var::rvx, side::right);
			double vmyl = ptrprb->side_val(SOL[cell], var::rvx, side::left);
			double vmyr = ptrprb->side_val(SOL[cell], var::rvx, side::right);
			double vrightl = ptrprb->side_val(SOL[cell+1], var::rvx, side::left);
			//double v49 = SOL[cell][0][var::rvx];

			double rleftr = ptrprb->side_val(SOL[cell - 1], var::r, side::right);
			double rmyl = ptrprb->side_val(SOL[cell], var::r, side::left);
			double rmyr = ptrprb->side_val(SOL[cell], var::r, side::right);
			double rrightl = ptrprb->side_val(SOL[cell + 1], var::r, side::left);
			
			//Взвешиваем скорости по Роу
			if (sqrt(rmyl)*vmyl + sqrt(rleftr)*vleftr > 0)
			{
				/*double ll = ptrprb->side_val(SOL[cell], sens, side::left);
				double rr = ptrprb->side_val(SOL[cell-1], sens, side::right);
				Ind[cell] += fabs(ll - rr);*/				
				Ind[cell] += fabs(ptrprb->side_val(SOL[cell], sens, side::left) - ptrprb->side_val(SOL[cell - 1], sens, side::right));
			}
			if (sqrt(rmyr)*vmyr + sqrt(rrightl)*vrightl < 0)
			{
				/*double ll = ptrprb->side_val(SOL[cell], sens, side::right);
				double rr = ptrprb->side_val(SOL[cell + 1], sens, side::left);
				Ind[cell] += fabs(ll - rr);*/
				Ind[cell] += fabs(ptrprb->side_val(SOL[cell], sens, side::right) - ptrprb->side_val(SOL[cell + 1], sens, side::left));
			}
			
			/*
			if (Ind[cell] == 0)
			{
				//double r48r = ptrprb->side_val(SOL[cell-1], sens, side::right);
				//double r49l = ptrprb->side_val(SOL[cell], sens, side::left);				
				//double r49r = ptrprb->side_val(SOL[cell], sens, side::right);
				//double r50l = ptrprb->side_val(SOL[cell + 1], sens, side::left);
				
				Ind[cell] += SOL[cell][0][var::rvx]>0 ? \
					fabs(ptrprb->side_val(SOL[cell], sens, side::left) - ptrprb->side_val(SOL[cell - 1], sens, side::right)) : \
					fabs(ptrprb->side_val(SOL[cell], sens, side::right) - ptrprb->side_val(SOL[cell + 1], sens, side::left));
			}
			*/
		}

		if (dynamic_cast<const ProblemMHD1D*>(ptrprb))
		{


			double vleftr = ptrprb->side_val(SOL[cell - 1], var::rvx, side::right);
			double vmyl = ptrprb->side_val(SOL[cell], var::rvx, side::left);
			double vmyr = ptrprb->side_val(SOL[cell], var::rvx, side::right);
			double vrightl = ptrprb->side_val(SOL[cell + 1], var::rvx, side::left);
			//double v49 = SOL[cell][0][var::rvx];

			double rleftr = ptrprb->side_val(SOL[cell - 1], var::r, side::right);
			double rmyl = ptrprb->side_val(SOL[cell], var::r, side::left);
			double rmyr = ptrprb->side_val(SOL[cell], var::r, side::right);
			double rrightl = ptrprb->side_val(SOL[cell + 1], var::r, side::left);

			//Взвешиваем скорости по Роу
			if (sqrt(rmyl)*vmyl + sqrt(rleftr)*vleftr > 0)
			{
				/*double ll = ptrprb->side_val(SOL[cell], sens, side::left);
				double rr = ptrprb->side_val(SOL[cell-1], sens, side::right);
				Ind[cell] += fabs(ll - rr);*/
				Ind[cell] += fabs(ptrprb->side_val(SOL[cell], sens, side::left) - ptrprb->side_val(SOL[cell - 1], sens, side::right));
			}
			if (sqrt(rmyr)*vmyr + sqrt(rrightl)*vrightl < 0)
			{
				/*double ll = ptrprb->side_val(SOL[cell], sens, side::right);
				double rr = ptrprb->side_val(SOL[cell + 1], sens, side::left);
				Ind[cell] += fabs(ll - rr);*/
				Ind[cell] += fabs(ptrprb->side_val(SOL[cell], sens, side::right) - ptrprb->side_val(SOL[cell + 1], sens, side::left));
			}

			/*
			if (Ind[cell] == 0)
			{
			//double r48r = ptrprb->side_val(SOL[cell-1], sens, side::right);
			//double r49l = ptrprb->side_val(SOL[cell], sens, side::left);
			//double r49r = ptrprb->side_val(SOL[cell], sens, side::right);
			//double r50l = ptrprb->side_val(SOL[cell + 1], sens, side::left);

			Ind[cell] += SOL[cell][0][var::rvx]>0 ? \
			fabs(ptrprb->side_val(SOL[cell], sens, side::left) - ptrprb->side_val(SOL[cell - 1], sens, side::right)) : \
			fabs(ptrprb->side_val(SOL[cell], sens, side::right) - ptrprb->side_val(SOL[cell + 1], sens, side::left));
			}
			*/
		}

		Ind[cell] = Ind[cell]/zn;
    }
    Ind[0] = 1.0;
    Ind[nx - 1] = 1.0;
}