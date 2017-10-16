#ifndef PARAMS_H_
#define PARAMS_H_

#include "defs.h"

//Тип данных "имя задачи"
enum CaseName { Sod,  BrioWu, Lax, Rarefaction, leftWoodward, rightWoodward, Shock, TestSmallTravelingWave, TestSmooth, \
	LeftTriangle};

#include "Boundary.h"
#include "BoundaryWall.h"
#include "BoundaryPeriodic.h"
#include "BoundarySoft.h"

#include "Indicator.h"
#include "IndicatorKrivodonova.h"  //Индикатор Криводоновой
//#include "IndicatorHarten.h"       //Индикатор Хартена
#include "IndicatorEverywhere.h"   //Индикатор, который срабатывает всегда (автоматически вызывается в МКР)
#include "IndicatorNowhere.h"      //Индикатор, который не срабатывает никогда

#include "Limiter.h"
#include "LimiterFinDiff.h"      //Метод конечных разностей
#include "LimiterWENO.h"         //Лимитер WENO
#include "LimiterWENO_S.h"         //Лимитер WENO_S
#include "LimiterHWENO_SC.h"         //Лимитер HWENO_SC
#include "LimiterHWENO_SC_Char.h"         //Лимитер HWENO_SC_Char
#include "LimiterHWENO.h"         //Лимитер WENO
//#include "LimiterWENOH.h"        //Лимитер WENOH

#include "Timestep.h"       
#include "TimestepEuler.h"   //Явный метод Эйлера
#include "TimestepRK2.h"     //Метод Рунге-Кутты 2-го порядка
#include "TimestepRK2TVD.h"  //Метод Рунге-Кутты 2-го порядка (схема TVD)
#include "TimestepRK3TVD.h"  //Метод Рунге-Кутты 3-го порядка (схема TVD)

#include "Flux.h"            
//#include "FluxKIR.h"          //Поток КИР
#include "FluxHLL.h"          //Поток HLL   
#include "FluxHLLC.h"         //Поток HLLC
//#include "FluxHLLE.h"
#include "FluxLaxFriedrichs.h"//Поток Лакса-Фридрихса
#include "FluxGodunovType.h"  //Поток типа Годунова (phiGodunov, phiCentralDiff, phivanLeer)

#include "Problem.h"           //Физическая постановка задачи
#include "ProblemGas1D.h"      //Одномерная газовая динамика, ур-я Эйлера
#include "ProblemMHD1D.h"      //Одномерная MГД
#include "ProblemTransfer1D.h" //Одномерное уравнение переноса

#include <vector>
#include <functional>
#include <memory>

using namespace std;

//Макрос для инициализации начальных условий
//(исключительно чтобы покороче писать)
#define InitF(name, expr) function<double(const double)> name = [=](const double x) -> double { return (expr); }

//Макрос для формирования вектора начальных условий
#define InitVector(name, ...) vector<std::function<double(const double)>> name = __VA_ARGS__


class Params
{
private:
	//функция, возвращающая тождественный нуль
	std::function<double(const double)> zero;

public:
	//Базовые параметры
	BaseParams basic;

	//Указатели на Постановку задачи, ГУ, индикатор, лимитер, интегратор по времени, поток
	Problem   *ptrprb;
	Boundary  *ptrbnd;
	Indicator *ptrind;
	Limiter   *ptrlim;
	Timestep  *ptrtst;
	Flux      *ptrflx;

	//Установка базовых параметров (исключая начальные условия)
	void SetBasicParams(double L, double T, int nx, double Co, int deltacnt)
	{
		basic.L = L;
		basic.T = T;
		basic.nx = nx;
		basic.Co = Co;
		basic.deltacnt = deltacnt;
		basic.h = L / nx;
		basic.tau = Co * L / nx;
	}


	//Два шаблона для установки типа задачи - без аргументов и с аргументами
	//template <typename T>
	//void SetProblem() { ptrprb = new T(basic); }

	template <typename T, typename... Arg>
	void SetProblem(Arg... arg) { ptrprb = new T(basic, arg...); }

	//Два шаблона для установки типа ГУ - без аргументов и с аргументами
	template <typename T>
	void SetBoundaryCondition() { ptrbnd = new T(basic, *ptrprb); }

	template <typename T, typename... Arg>
	void SetBoundaryCondition(Arg... arg) { ptrbnd = new T(basic, *ptrprb, arg...); }

	//Два шаблона для установки типа индикатора - без аргументов и с аргументами
	template <typename T>
	void SetIndicator()	{ ptrind = new T(basic, *ptrprb); }

	template <typename T, typename... Arg>
	void SetIndicator(Arg... arg)
	{
		ptrind = new T(basic, *ptrprb, arg...);
	}

	//Два шаблона для установки типа лимитера - без аргументов и с аргументами
	template <typename T>
	void SetLimiter() { ptrlim = new T(basic, *ptrprb, *ptrind); }

	template <typename T, typename... Arg> 
	void SetLimiter(Arg... arg) { ptrlim = new T(basic, *ptrprb, *ptrind, arg...); }

	//Два шаблона для установки типа интегратора по времени - без аргументов и с аргументами
	template <typename T>
	void SetTimestep() { ptrtst = new T(basic, ptrprb->dim, ptrprb->nshape, *ptrbnd); }

	template <typename T, typename... Arg>
	void SetTimestep(Arg... arg) { ptrtst = new T(basic, ptrprb->dim, ptrprb->nshape, *ptrbnd, arg...); }

	//Два шаблона для установки типа численного потока - без аргументов и с аргументами
	template <typename T>
	void SetFlux() { ptrflx = new T(basic, *ptrprb); }		

	template <typename T, typename... Arg>
	void SetFlux(Arg... arg) { ptrflx = new T(basic, *ptrprb, arg...); }


	//Конструктор:
	//problem - идентификатор решаемой задачи
	Params(const CaseName SovingCase);

	//Деструктор
	~Params()
	{
		delete ptrind;
		delete ptrlim;
		delete ptrbnd;
		delete ptrflx;
		delete ptrtst;
	};

};

#endif
