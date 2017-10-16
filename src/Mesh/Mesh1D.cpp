#include <iostream>
#include <functional>

#include "Mesh1D.h"
#include "Integrator.h"


Mesh1D::Mesh1D(const Params& prm)
{
	ptrprm = &prm;

	int nx = ptrprm->basic.nx;
	int dim = ptrprm->ptrprb->dim;
	int nshape = ptrprm->ptrprb->nshape;

	SOL1.resize(nx + 2);
	SOL2.resize(nx + 2);

	for (int cell = 0; cell < (nx + 2); ++cell)
	{
		SOL1[cell].resize(nshape);
		SOL2[cell].resize(nshape);

		for (int shape = 0; shape < nshape; ++shape)
		{
			SOL1[cell][shape].resize(dim);
			SOL2[cell][shape].resize(dim);
		}
	}

    ptrSOL = &SOL2; ptrSOLnew = &SOL1;


    //Начальное условие - средние значения и наклоны
	step_num = 0;
	set_initial();
	
    telefile.open("../RKDG1D/telefile.txt");
    init_telefile();

	//Сохранение в файл начального состояния
	add_to_telefile();

    //telefile.close(); 
    //telefile.clear();

};

Mesh1D::~Mesh1D()
{
    telefile.close();
};


//Запонение вектора решения начальными условиями
void Mesh1D::set_initial()
{	
	//Для удобства ссылки:
	const double& h = ptrprm->basic.h;
	const int& nx = ptrprm->basic.nx;
	const int& dim = ptrprm->ptrprb->dim;
	const int nshape = ptrprm->ptrprb->nshape;

	
	int ngp = 3;
	Integrator GP(IntegrPoints::Gauss, ngp);
	
	auto& shapes = ptrprm->ptrprb->shapefunc;
	auto& shapenrm2 = ptrprm->ptrprb->shapefuncnorm2;

	for (int cell = 0; cell < nx; ++cell)
	{
		function<vector<double>(double)> initial = [&](double x) {return ptrprm->ptrprb->initial_var(cnt(cell) + 0.5*h*x); };
				
		for (int shape = 0; shape < nshape; ++shape)
			(*ptrSOLnew)[cell][shape] = GP.integrate([&](double pt) {return (initial(pt) * shapes[shape](pt)) / shapenrm2[shape]; }, dim);
	}

	/*
	for (int cell = 0; cell < nx; ++cell)	
	{
		function<vector<double>(double)> initial = [&](double x) {return ptrprm->ptrprb->initial_var(cnt(cell) + 0.5*h*x);};
			
		vector<double> intf0 = GP.integrate([&](double pt) {return initial(pt);}, dim);		
		
		(*ptrSOLnew)[cell][0] = intf0;
		
		if (nshape > 1)
		{
			vector<double> intf1 = GP.integrate([&](double pt) {return initial(pt) * pt * 3.0; }, dim);
			(*ptrSOLnew)[cell][1] = intf1;
		}

		if (nshape > 2)
		{
			vector<double> intf2 = GP.integrate([&](double pt) {return initial(pt) * (1.5 * sqr(pt) - 0.5) * 5.0; }, dim);
			(*ptrSOLnew)[cell][2] = intf2;
		}

	}
	*/
}

//Выполнение одного шага по времени
void Mesh1D::TimeStepping()
{
	//Для удобства ссылки:
	const double& tau = ptrprm->basic.tau;
	const int& deltacnt = ptrprm->basic.deltacnt;
	
	//Инкремент номера шага по времени
	step_num++;

	//Техническая процедура переприсваивания указателей
	//Обмен указателей на решения/наклоны на предыдущем и текущем шагах,
	//выполняется в начале каждого шага
	std::swap(ptrSOL, ptrSOLnew);

	//Выполнение интегрирования на одном шаге выбранным численным методом
	ptrprm->ptrtst->runstep(SOL(), SOLnew(), tau, *ptrprm->ptrflx, *ptrprm->ptrlim);
	
	/*
	if (step_num==5000)
	for (int i = 0; i < ptrprm->basic.nx; i++)
	{
		telefile << Unew()[i][0] << endl;
	}
	//*/


	//Сохранение в файл, если номер шага делится нацело на deltacnt
	if (!(step_num%deltacnt))
	{
		add_to_telefile();
		cout << "t = " << step_num * tau << endl;
	}


};



//Инициализация телефайла
void Mesh1D::init_telefile()
{
	//Для удобства ссылки:
	//*
	const double& h = ptrprm->basic.h;
	const double& tau = ptrprm->basic.tau;
	const int& nx = ptrprm->basic.nx;
	const int& deltacnt = ptrprm->basic.deltacnt;
	const int& nshape = ptrprm->ptrprb->nshape;

	telefile.precision(13);
	telefile << "Nx= " << nx << " h= " << h << " tau= " << tau << " deltacnt= " << deltacnt << " nshape= "<< nshape << endl ;
	for (int cell = 0; cell < nx; ++cell)
        telefile << cnt(cell) << " ";
    telefile << endl;
	//*/
}

//Сохранение телеметрии
void Mesh1D::add_to_telefile()
{	
	//*
	const vector<vector<vector<double>>>& SOL = SOLnew();

	const int& nx = ptrprm->basic.nx;
	const int& dim = ptrprm->ptrprb->dim;
	const int nshape = ptrprm->ptrprb->nshape;

	for (int cell = 0; cell < nx; ++cell)
	{
		for (var val = (var)0; val < dim; val = (var)(val + 1))
		{
			for (int shape = 0; shape < nshape; ++shape)
			{
				telefile << SOL[cell][shape][val] << " ";
			}
		}

		ptrprm->ptrprb->more_to_file(telefile, SOL, cell);
	}

	

    telefile << endl;
	//*/

}//add_to_telefile;

