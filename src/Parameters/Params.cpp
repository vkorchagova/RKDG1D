#include "Params.h"
#include "defs.h"

Params::Params(const CaseName SovingCase)
{	
	zero = [](const double) -> double { return 0.0; };
	
	double L; double T; int nx; double Co; int deltacnt;
	int nshape;
	
	//Заполнение свойств объекта Params исходя из математической постановки задачи
	//и параметров дискретизации
	switch (SovingCase)
	{

	case Sod:
	{
		L = 1.0;
        T = 0.2;

        nx = 200;
		Co = 0.1;
		deltacnt = 1;
		SetBasicParams(L, T, nx, Co, deltacnt);		//Обязательная команда

		
		//Специфические параметры задачи
		double gamma = 1.4; //показатель адиабаты

		//Задание начальных условий
		InitF(rho, (x <= 0.5*L) ? 1.0 : 0.125);
		InitF(p, (x <= 0.5*L) ? 1.0 : 0.1);
		InitVector(initv, { rho, zero, zero, zero, p });
		
		//Задание решаемой задачи, ГУ, индикатора, лимитера, схемы интегрирования, потока
		//Тип - в угловых скобках, параметры (произвольное число, понятное конструктору) - в круглых
        nshape = 3;
		SetProblem<Gas1D>(nshape, gamma, initv);
		SetBoundaryCondition<Wall>();
		//SetFlux<HLLC>(RoeSoundVel);
        SetFlux<LaxFriedrichs>(RoeSoundVel);
		//SetFlux<KIR>();
		SetIndicator<Krivodonova>(var::r);
        //SetIndicator<Harten>(var::r, 2.0);
        //SetIndicator<Nowhere>();
        SetLimiter<WENO>(2.0);            //Инициализируется после индикатора
		//SetLimiter<FinDiff>();
        SetTimestep<RK3TVD>();			  //Инициализируется после гран.условий
	}
	break;

	case BrioWu:
	{
		L = 1.0;
		T = 0.1;

		nx = 100;
		Co = 0.1;
		deltacnt = 1;
		SetBasicParams(L, T, nx, Co, deltacnt);		//Обязательная команда


		//Специфические параметры задачи
		double gamma = 2.0; //показатель адиабаты
		double Hx = 0.75;

		//Задание начальных условий
		InitF(rho, (x <= 0.5*L) ? 1.0 : 0.125);
		InitF(p, (x <= 0.5*L) ? 1.0 : 0.1);
		InitF(Hy, (x <= 0.5*L) ? 1.0 : -1.0);
		InitVector(initv, { rho, zero, zero, zero, p, Hy, zero });

		//Задание решаемой задачи, ГУ, индикатора, лимитера, схемы интегрирования, потока
		//Тип - в угловых скобках, параметры (произвольное число, понятное конструктору) - в круглых
        nshape = 1;
		SetProblem<MHD1D>(nshape, gamma, Hx, initv);
		SetBoundaryCondition<Wall>();
		//SetFlux<HLLC>(RoeSoundVel);
		SetFlux<LaxFriedrichs>(RoeSoundVel);
		//SetFlux<KIR>();
		SetIndicator<Krivodonova>(var::r);
		//SetIndicator<Harten>(var::r, 2.0);
		//SetIndicator<Nowhere>();
		SetLimiter<WENO_S>(2.0);            //Инициализируется после индикатора
		//SetLimiter<FinDiff>();
		SetTimestep<RK2TVD>();			  //Инициализируется после гран.условий
	}
	break;

	case Lax: // (см. статью по HWENO_SC)
	{
		L = 10.0;
		T = 1.3;

		nx = 100;
		Co = 0.1;
		deltacnt = 10;
		SetBasicParams(L, T, nx, Co, deltacnt);		//Обязательная команда


		//Специфические параметры задачи
		double gamma = 1.4; //показатель адиабаты

		//Задание начальных условий
		InitF(rho, (x <= 0.5*L) ? 0.445 : 0.5);
		InitF(p, (x <= 0.5*L) ? 3.528 : 0.571);
		InitF(u, (x <= 0.5*L) ? 0.698 : 0.0);
		InitVector(initv, { rho, u, zero, zero, p });

		//Задание решаемой задачи, ГУ, индикатора, лимитера, схемы интегрирования, потока
		//Тип - в угловых скобках, параметры (произвольное число, понятное конструктору) - в круглых
		nshape = 3;
		SetProblem<Gas1D>(nshape, gamma, initv);
		SetBoundaryCondition<Soft>();
		SetFlux<HLLC>(RoeSoundVel);
		//SetFlux<LaxFriedrichs>(RoeSoundVel);
		//SetFlux<KIR>();
		SetIndicator<Krivodonova>(var::r);
		//SetIndicator<Harten>(var::r, 2.0);
		//SetIndicator<Nowhere>();
		SetLimiter<WENO>(2.0);            //Инициализируется после индикатора
		//SetLimiter<FinDiff>();
		SetTimestep<RK3TVD>();			  //Инициализируется после гран.условий
	}
	break;


	//case Rarefaction:
	//{		
	//	L = 1.0;
	//	T = 0.005;

	//	nx = 100;
	//	Co = 0.001;
	//	deltacnt = 100;
	//	SetBasicParams(L, T, nx, Co, deltacnt);		//Обязательная команда

	//	//Специфические параметры задачи
	//	double gamma = 1.4; //показатель адиабаты

	//	//Задание начальных условий
	//	InitF(rho, 1.0);
	//	InitF(u, (x < 0.5*L) ? -2.0 : 2.0);
	//	InitF(p, 0.4);
	//	InitVector(initv, { rho, u, zero, zero, p });

	//	//Задание решаемой задачи, ГУ, индикатора, лимитера, схемы интегрирования, потока
	//	//Тип - в угловых скобках, параметры (произвольное число, понятное конструктору) - в круглых
	//	nshape = 3;
	//	SetProblem<Gas1D>(nshape, gamma, initv);
	//	SetBoundaryCondition<Soft>();
	//	SetFlux<HLLC>(RoeSoundVel);
	//	//SetFlux<LaxFriedrichs>(RoeSoundVel);
	//	//SetFlux<KIR>();
	//	SetIndicator<Krivodonova>(var::r);
	//	//SetIndicator<Harten>(var::r, 2.0);
	//	//SetIndicator<Nowhere>();
	//	SetLimiter<WENO>(2.0);            //Инициализируется после индикатора
	//	//SetLimiter<FinDiff>();
	//	SetTimestep<RK3TVD>();			  //Инициализируется после гран.условий
	//}
	//break;



    case leftWoodward:
    {
        L = 1.0;
        T = 0.012;
            
        nx = 800;
        Co = 0.0001;
        deltacnt = 100;
        SetBasicParams(L, T, nx, Co, deltacnt); 		//Обязательная команда
            
        //Специфические параметры задачи
        double gamma = 1.4; //показатель адиабаты
            
        //Задание начальных условий
        InitF(rho, 1.0);
        InitF(p, (x <= 0.5*L) ? 1000.0 : 0.01);
        InitVector(initv, { rho, zero, zero, zero, p });
            
        //Задание решаемой задачи, ГУ, индикатора, лимитера, схемы интегрирования, потока
        //Тип - в угловых скобках, параметры (произвольное число, понятное конструктору) - в круглых
		nshape = 2;
		SetProblem<Gas1D>(nshape, gamma, initv);
		SetBoundaryCondition<Wall>();
		SetFlux<HLLC>(RoeSoundVel);
		//SetFlux<LaxFriedrichs>(RoeSoundVel);
		//SetFlux<KIR>();
		SetIndicator<Krivodonova>(var::r);
		//SetIndicator<Harten>(var::r, 2.0);
		//SetIndicator<Everywhere>();
		SetLimiter<WENO>(2.0);            //Инициализируется после индикатора
		//SetLimiter<FinDiff>();
		SetTimestep<RK3TVD>();			  //Инициализируется после гран.условий
    }
    break;
          



    case rightWoodward:
    {
        L = 1.0;
        T = 0.01;
            
        nx = 200;
        Co = 0.01;
        deltacnt = 1;
        SetBasicParams(L, T, nx, Co, deltacnt); 		//Обязательная команда
            
        //Специфические параметры задачи
        double gamma = 1.4; //показатель адиабаты
            
        //Задание начальных условий
		//InitF(rho, (x <= 0.5*L) ? 1.0 : 8.0);
        //InitF(p, (x <= 0.5*L) ? 1.0 : 10.0);


		InitF(rho, (x <= 0.5*L) ? 1.0 : 1.0);
		InitF(p, (x <= 0.5*L) ? 0.01 : 100.0);

        InitVector(initv, { rho, zero, zero, zero, p });
            
        //Задание решаемой задачи, ГУ, индикатора, лимитера, схемы интегрирования, потока
        //Тип - в угловых скобках, параметры (произвольное число, понятное конструктору) - в круглых
		nshape = 3;
		SetProblem<Gas1D>(nshape, gamma, initv);
		SetBoundaryCondition<Soft>();
		SetFlux<HLLC>(RoeSoundVel);
		//SetFlux<LaxFriedrichs>(RoeSoundVel);
		//SetFlux<KIR>();
		SetIndicator<Krivodonova>(var::r);
		//SetIndicator<Harten>(var::r, 2.0);
		//SetIndicator<Nowhere>();
		SetLimiter<HWENO_SC>(2.0);            //Инициализируется после индикатора
		//SetLimiter<FinDiff>();
		SetTimestep<RK3TVD>();			  //Инициализируется после гран.условий
    }
    break;




    case Shock:
    {
        L = 1.0;
        T = 0.03;
            
        nx = 100;
        Co = 0.01;
        deltacnt = 1;
        SetBasicParams(L, T, nx, Co, deltacnt); 		//Обязательная команда
            
        //Специфические параметры задачи
        double gamma = 1.4; //показатель адиабаты
            
        //Задание начальных условий
        InitF(rho, (x <= 0.5*L) ? 5.99924 : 5.99242);
        InitF(u, (x < 0.5*L) ? 19.5975 : -6.19633);
        InitF(p, (x <= 0.5*L) ? 460.894 : 46.0950);
        InitVector(initv, { rho, u, zero, zero, p });
            
        //Задание решаемой задачи, ГУ, индикатора, лимитера, схемы интегрирования, потока
        //Тип - в угловых скобках, параметры (произвольное число, понятное конструктору) - в круглых
		nshape = 3;
		SetProblem<Gas1D>(nshape, gamma, initv);
		SetBoundaryCondition<Soft>();
		SetFlux<HLLC>(RoeSoundVel);
		//SetFlux<LaxFriedrichs>(RoeSoundVel);
		//SetFlux<KIR>();
		SetIndicator<Krivodonova>(var::r);
		//SetIndicator<Harten>(var::r, 2.0);
		//SetIndicator<Nowhere>();
		SetLimiter<HWENO_SC>(2.0);            //Инициализируется после индикатора
		//SetLimiter<FinDiff>();
		SetTimestep<RK3TVD>();			  //Инициализируется после гран.условий
    }
    break;


    case TestSmallTravelingWave:
    {
		const double A = 1e-6;
		const double PI = 3.1415926535897932384626433832795;

		L = 1.0;
		T = 10.0;

		nx = 10;
		Co = 0.1;
		deltacnt = 1;// nx * 10;
		SetBasicParams(L, T, nx, Co, deltacnt); 		//Обязательная команда

		//Специфические параметры задачи
		double gamma = 1.4; //показатель адиабаты

		//Задание начальных условий
		InitF(rho, 1.0 + A * sin(2.0 * PI *x));		
		InitF(p, 1.0 / gamma);
		InitF(u, -1.0);
		InitVector(initv, { rho, u, zero, zero, p });
		
		//Задание решаемой задачи, ГУ, индикатора, лимитера, схемы интегрирования, потока
		//Тип - в угловых скобках, параметры (произвольное число, понятное конструктору) - в круглых
		nshape = 3;
		SetProblem<Gas1D>(nshape, gamma, initv);
		SetBoundaryCondition<Periodic>();
		SetFlux<HLLC>(RoeSoundVel);
		
		//SetIndicator<Krivodonova>(var::r);
		//SetLimiter<FinDiff>();           //Инициализируется после индикатора
		
		SetIndicator<Nowhere>();
		//SetIndicator<Harten>(var::r, 2.0);
		SetLimiter<FinDiff>();             //Инициализируется после индикатора
		SetTimestep<RK3TVD>();			   //Инициализируется после гран.условий
	}
	break;

/*
	case TestSmooth:
	{
		const double A = 500.0;

		L = 1.0;
		T = 0.2;

		nx = 100;
		Co = 0.1;
		deltacnt = 10;
		SetBasicParams(L, T, nx, Co, deltacnt); 		//Обязательная команда

		//Специфические параметры задачи
		double gamma = 1.4; //показатель адиабаты

		//Задание начальных условий
		InitF(rho, ((x >= 0.25*L) && (x <= 0.75*L)) ? 1.0 + A * sqr((x - 0.25*L)*(x - 0.75*L)) : 1.0);
		InitF(drho, ((x >= 0.25*L) && (x <= 0.75*L)) ? 2.0 * A * (x - 0.25*L)*(x - 0.75*L)*(2.0*x - 1.0) : 0.0);
		InitF(p, ((x >= 0.25*L) && (x <= 0.75*L)) ? (1.0 + A * sqr((x - 0.25*L)*(x - 0.75*L)))*0.4 : 1.0 * 0.4);
		InitF(dp, ((x >= 0.25*L) && (x <= 0.75*L)) ? (2.0 * A * (x - 0.25*L)*(x - 0.75*L)*(2.0 * x - 1.0))*0.4 : 0.0*0.4);
		InitVector(initv, { rho, zero, zero, zero, p, drho, zero, zero, zero, dp });

		//Задание решаемой задачи, ГУ, индикатора, лимитера, схемы интегрирования, потока
		//Тип - в угловых скобках, параметры (произвольное число, понятное конструктору) - в круглых
		SetProblem<Gas1D>(gamma, initv);
		SetBoundaryCondition<Wall>();
		SetFlux<HLLC>(RoeSoundVel);                   
		SetIndicator<Nowhere>();
		SetLimiter<WENOH>(1.5);             //Инициализируется после индикатора
		SetTimestep<RK2>();					//Инициализируется после гран.условий	
	}
	break;
*/

/*
	case LeftTriangle:
	{
		L = 1.0;
		T = 0.5;

		nx = 200;
		Co = 0.1;
		deltacnt = 1;
		
		SetBasicParams(L, T, nx, Co, deltacnt); 		//Обязательная команда

		//Задание начальных условий
		InitF(U, ((x <= 0.50*L) && (x >= 0.25*L)) ? (x - 0.25) / 0.25 : 0.0);
		InitF(dU, ((x <= 0.50*L) && (x >= 0.25*L)) ? 4.0 : 0.0);
		InitVector(initv, { U, dU });

		//Задание решаемой задачи, ГУ, индикатора, лимитера, схемы интегрирования, потока
		//Тип - в угловых скобках, параметры (произвольное число, понятное конструктору) - в круглых
		
		//Линейное уравнение переноса
		//SetProblem<Transfer1D>(initv, [](const double U){ return U; });
		
		// Квазилинейное уравнение перенос
		SetProblem<Transfer1D>(initv, [](const double U){ return 0.5*U*U; });
		
		SetBoundaryCondition<Wall>();
		SetFlux<GodunovType>(vanLeer, 1e-10);  
		SetIndicator<Harten>(var::U, 2.0);
		SetLimiter<WENO>(2.0);                  //Инициализируется после индикатора
		SetTimestep<RK2>();						//Инициализируется после гран.условий
	}
	break;
*/

	}//switch
	
}

