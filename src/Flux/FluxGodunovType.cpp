#include "FluxGodunovType.h"



//Нелинейные функции-ограничители при расчете численного потока:
//  - поток С.К. Годунова (против потока, монотонна, 1-й порядок)
//  - поток, рассчитанный по полусумме решений в левой и правой ячейках (немонотонна, 2-й порядок)
//  - поток ван-Лира (монотонна, 2-й порядок)

double phiGodunov(const double r) { return 0.0; };
double phiCentralDiff(const double r) { return 1.0; }
double phivanLeer(const double r) { return (r + fabs(r)) / (1.0 + fabs(r)); };


//КОНСТРУКТОР Годунова
FluxGodunovType::FluxGodunovType(const BaseParams& prm, Problem& prb,\
	double(*const phi_)(const double), double thres_)
	: Flux(prm, prb), phi(phi_), thres(thres_)
{
}//FluxGodunovType::FluxGodunovType

//ДЕСТРУКТОР Годунова
FluxGodunovType::~FluxGodunovType()
{}; 

double FluxGodunovType::FluxNum(const vector<double>& fluxL, const vector<double>& fluxR, const double r, const var val)
{
	const vector<double>& upwind = fluxL;

	double p = phi(r);

	return (1.0 - p)*upwind[val] + p*0.5*(fluxL[val]+fluxR[val]);
}


void FluxGodunovType::step(const vector<vector<vector<double>>>& SOL, \
	vector<vector<vector<double>>>& DSOL, \
	const double cft)
{
	int nx = ptrprm->nx;
	double h = ptrprm->h;

	//Находим конвективные потоки
	ptrprb->convFlux(SOL);

	double rleft = 0.0, rright = 0.0;
	double flL = 0.0, flR = 0.0;
	vector<double> DSOLpred (ptrprb->nshape, 0.0);	

	//Основной цикл по ячейкам
	for (int cell = 0; cell < nx; ++cell)
	{
		//Установка ссылок на решения на трех ячейках	:
		//  myu(), myv() - на своей ячейке
		//  leftu(), leftv() - на соседней слева ячейке
		//  rightu(), rightv() - на соседней справа ячейке 
		//и на конв.потоки на трех ячейках:
		//  myflux(), myfluxL(), myfluxR() - конв.потоки на своей ячейке (по центру, слева, справа)
		//  rightfluxL(), leftfluxR() - конв.потоки на правой и левой соседних ячейках (слева, справа)
		setlocsolflux(SOL, cell);

		//Пересчитываем средние значения и потоки
		for (int val = 0; val < ptrprb->dim; ++val)
		{
			if (1>0)
			{
				rleft = (side_val(leftsol(),     (var)val, side::left) - \
					     side_val(leftleftsol(), (var)val, side::right)) / \
					    (side_val(mysol(),       (var)val, side::left) - \
					     side_val(leftleftsol(), (var)val, side::right) + 1e-10);

				rright = (side_val(mysol(),    (var)val, side::left) - \
					      side_val(leftsol(),  (var)val, side::right)) / \
					     (side_val(rightsol(), (var)val, side::left) - \
			  		      side_val(mysol(),    (var)val, side::right) + 1e-10);
			}
			
			//ТРЕБУЕТ ПЕРЕРАБОТКИ!!!
			if (1<0)
			{				
				cout << "NOT_MODIFIED!!!" << endl;
				exit(1);		
			}

			flR = FluxNum(myfluxR(),   rightfluxL(), rright, (var)val);
			flL = FluxNum(leftfluxR(), myfluxL(),    rleft,  (var)val);

			DSOL[cell][0][val] = -(cft / h)*(flR - flL);

            //По формуле центральных прямоугольников
			//DVpred = -3.0*(cft / h)*(-2.0*(myflux[val]) + flR + flL);
            
            //По формуле Симпсона
			if (ptrprb->nshape > 1)
			{
				DSOLpred[1] = -3.0*(cft / h)*(-2.0*(myfluxL()[val] + 4.0*myflux()[val] + myfluxR()[val]) / 6.0 + flR + flL);
				DSOL[cell][1][val] = (fabs(DSOLpred[1]) < thres) ? 0.0 : DSOLpred[1];
			}
			
			if (ptrprb->nshape > 2)
			{
				DSOLpred[2] = -5.0*(cft / h)*(-1.0*(-myfluxL()[val] + myfluxR()[val]) + flR - flL);
				DSOL[cell][2][val] = (fabs(DSOLpred[2]) < thres) ? 0.0 : DSOLpred[2];
			}			
		}

	}//for cell

}//FluxGodunovType::step

void FluxGodunovType::setlocsolflux(const vector<vector<vector<double>>>& SOL, const int cell)
{
	Flux::setlocsolflux(SOL, cell);
	
	int nx = ptrprm->nx;

	ptr_leftleftsol = (cell > 1) ? &(SOL[cell - 2]) : &(SOL[nx]);	
}

