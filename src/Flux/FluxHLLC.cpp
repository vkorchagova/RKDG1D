#include "FluxHLLC.h"
#include "Integrator.h"

//КОНСТ� УКТО�  HLLC
FluxHLLC::FluxHLLC(const BaseParams& prm, Problem& prb, SoundVelType soundvel) : Flux(prm, prb), SoundVel(soundvel)
{
    int nx = ptrprm->nx;

    //Заготовки под L
    L.resize(nx + 1);
    for (int row = 0; row < nx + 1; ++row)
    {
        L[row].resize(2);
    }

    //Заготовки под необходимые (временные) векторы
    LdUast.resize(5), RdUast.resize(5), hllcL.resize(5), hllcR.resize(5), UastL.resize(5), UastR.resize(5);
}

//ДЕСТ� УКТО�  HLLC
FluxHLLC::~FluxHLLC()
{};


void FluxHLLC::getUast(const int cell, const double Sast, \
	const vector<vector<double>>& leftsol, \
	const vector<vector<double>>& mysol, \
	vector<double>& Uast)
{
    if (Sast >= 0)
    {
		double mnog = side_val(leftsol, var::r, side::right) * \
			(L[cell][0] - side_val(leftsol, var::vx, side::right)) / (L[cell][0] - Sast);
        
        Uast[0] = 1.0;
        Uast[1] = Sast;
		Uast[2] = side_val(leftsol, var::vy, side::right);
		Uast[3] = side_val(leftsol, var::vz, side::right);
		Uast[4] = side_val(leftsol, var::e, side::right) / side_val(leftsol, var::r, side::right) + \
			(Sast - side_val(leftsol, var::vx, side::right))* \
			(Sast + side_val(leftsol, var::p, side::right) / \
			(side_val(leftsol, var::r, side::right)*(L[cell][0] - side_val(leftsol, var::vx, side::right))));

        for (int val = 0; val < 5; ++val)
            Uast[val] *= mnog;
    }
    else //if Sast < 0
    {
        double mnog = side_val(mysol, var::r, side::left) * \
			(L[cell][1] - side_val(mysol, var::vx, side::left)) / (L[cell][1] - Sast);

        Uast[0] = 1.0;
        Uast[1] = Sast;
		Uast[2] = side_val(mysol, var::vy, side::left);
		Uast[3] = side_val(mysol, var::vz, side::left);
		Uast[4] = side_val(mysol, var::e, side::left) / side_val(mysol, var::r, side::left) + \
			(Sast - side_val(mysol, var::vx, side::left))* \
			(Sast + side_val(mysol, var::p, side::left) / \
			(side_val(mysol, var::r, side::left)*(L[cell][1] - side_val(mysol, var::vx, side::left))));
        
        for (int val = 0; val < 5; ++val)
            Uast[val] *= mnog;
    }//end else
}


double FluxHLLC::getSast(const int cell, \
	const vector<vector<double>>& leftsol, \
	const vector<vector<double>>& mysol)
{
	return (side_val(mysol, var::p, side::left) - side_val(leftsol, var::p, side::right) + \
		side_val(leftsol, var::rvx, side::right) * (L[cell][0] - side_val(leftsol, var::vx, side::right)) - \
		side_val(mysol, var::rvx, side::left) * (L[cell][1] - side_val(mysol, var::vx, side::left))) / \
		(side_val(leftsol, var::r, side::right) * (L[cell][0] - side_val(leftsol, var::vx, side::right)) - \
		side_val(mysol, var::r, side::left) * (L[cell][1] - side_val(mysol, var::vx, side::left)));
}


void FluxHLLC::step(const vector<vector<vector<double>>>& SOL, \
	vector<vector<vector<double>>>& DSOL, \
	const double cft)
{
    int nx = ptrprm->nx;
    double h = ptrprm->h;
	int dim = ptrprb->dim;

	//Находим конвективные потоки
	ptrprb->convFlux(SOL);

	

    //Находим собственные числа
	ptrprb_toGas->lambda(SOL, SoundVel, L, { 0, dim - 1 });
    
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

		/*
		if (fabs(L[cell][1]) < 1e-10)
		L[cell][1] = 0.0;

		if (fabs(L[cell+1][1]) < 1e-10)
		L[cell+1][1] = 0.0;
		*/

		/*
		if ((L[cell][1]) < 0)
		{
		cout << cell << " " << L[cell][1] << endl;
		cin.get();
		}
		*/

		double SastL = getSast(cell, leftsol(), mysol());
		double SastR = getSast(cell + 1, mysol(), rightsol());

		getUast(cell, SastL, leftsol(), mysol(), UastL);
		getUast(cell + 1, SastR, mysol(), rightsol(), UastR);

		//Проверка направлений переноса
		if (L[cell][0] > 0)
		for (size_t val = 0; val < dim; ++val)
			hllcL[val] = leftfluxR()[val];
		else
		if (SastL >= 0)
		for (size_t val = 0; val < dim; ++val)
		{
			LdUast[val] = UastL[val] - side_val(leftsol(), (var)val, side::right);
			hllcL[val] = leftfluxR()[val] + L[cell][0] * LdUast[val];
		}
		else
		if (L[cell][1] > 0)
		for (size_t val = 0; val < dim; ++val)
		{
			LdUast[val] = UastL[val] - side_val(mysol(), (var)val, side::left);
			hllcL[val] = myfluxL()[val] + L[cell][1] * LdUast[val];
		}
		else
		for (size_t val = 0; val < dim; ++val)
			hllcL[val] = myfluxL()[val];


		//Проверка направлений переноса
		if (L[cell + 1][0] > 0)
		for (size_t val = 0; val < dim; ++val)
			hllcR[val] = myfluxR()[val];
		else
		if (SastR >= 0)
		for (size_t val = 0; val < dim; ++val)
		{
			RdUast[val] = UastR[val] - side_val(mysol(), (var)val, side::right);
			hllcR[val] = myfluxR()[val] + L[cell + 1][0] * RdUast[val];
		}
		else
		if (L[cell + 1][1] > 0)
		for (size_t val = 0; val < dim; ++val)
		{
			RdUast[val] = UastR[val] - side_val(rightsol(), (var)val, side::left);
			hllcR[val] = rightfluxL()[val] + L[cell + 1][1] * RdUast[val];
		}
		else
		for (size_t val = 0; val < dim; ++val)
			hllcR[val] = rightfluxL()[val];


		int ngp = 3;
		Integrator GP(IntegrPoints::Gauss, ngp);
		
		function<vector<double>(double)> flx = [&](double x) {return ptrprb->getFlux(ptrprb->gauss_val(mysol(), x));};
		
        vector<double> intFL1 = GP.integrate([&](double pts) {return flx(pts)*2.0;}, dim);
        vector<double> intFL2 = GP.integrate([&](double pts) {return flx(pts)*6.0*pts;}, dim);0;
				
        //Пересчитываем средние значения и потоки
        for (size_t val = 0; val < 5; ++val)
        {
            DSOL[cell][0][val] = -(cft / h) * (hllcR[val] - hllcL[val]);
            DSOL[cell][1][val] = -3.0*(cft / h) * (-intFL1[val] + hllcR[val] + hllcL[val]);		

			if (ptrprb->nshape > 2)
				DSOL[cell][2][val] = -5.0*(cft / h) * (-intFL2[val] + hllcR[val] - hllcL[val]);
        }
    }
}
