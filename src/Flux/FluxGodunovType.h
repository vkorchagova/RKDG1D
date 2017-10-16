#ifndef FLUXGODUNOVTYPE_H_
#define FLUXGODUNOVTYPE_H_

#include <vector>
#include <cmath>
#include <algorithm>

#include "defs.h"
#include "Flux.h"
#include "Problem.h"
#include "ProblemTransfer1D.h"


using namespace std;

#define GodunovType FluxGodunovType


//Ќелинейные функции-ограничители при расчете численного потока:
//  - поток —. . √одунова (против потока, монотонна, 1-й пор€док)
//  - поток, рассчитанный по полусумме решений в левой и правой €чейках (немонотонна, 2-й пор€док)
//  - поток ван-Ћира (монотонна, 2-й пор€док)

#define Godunov phiGodunov
double phiGodunov(const double r);

#define CentralDiff phiCentralDiff
double phiCentralDiff(const double r);

#define vanLeer phivanLeer
double phivanLeer(const double r);


class FluxGodunovType :
	public Flux
{
private:

	//ћалое число - отсечка дл€ наклонов:
	double thres;

	//—сылка на функцию-ограничитель (phiGodunov, phivanLeer...)
	double(*const phi)(const double);

	//‘ункци€ численного потока, 
	//вычисл€емого по решению слева (UL) и справа (UR) от разрыва,
	//значению r дл€ переменной val
	double FluxNum(const vector<double>& fluxL, const vector<double>& fluxR, const double r, const var val);

	//”казатели на векторы решени€ в €чейке через одну слева
	const vector<vector<double>> *ptr_leftleftsol;

protected:
	const vector<vector<double>>& leftleftsol()  { return *ptr_leftleftsol; };

	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я
	//”становка приватных указателей (чтобы работали вышеперечисленные ссылки): 
	// - на конв. потоки в своей €чейке
	// - на предельные (слева и справа) потоки в фиктивных €чейках
	// - на решение в своей €чейке и в соседних €чейках
	virtual void setlocsolflux(const vector<vector<vector<double>>>& SOL, const int cell);


public:
	// онструктор (prm - объект, содержащий параметры задачи, 
	//            phi - функци€-ограничитель)
	FluxGodunovType(const BaseParams& prm, Problem& prb,\
		double(*const phi_)(const double), double thres);

	//ƒеструктор
	~FluxGodunovType();

	//–≈јЋ»«ј÷»я ¬»–“”јЋ№Ќќ… ‘”Ќ ÷»» 
	//—обственно, шаг расчета
	//заполнение приращений решений (DU) и наклонов (DV) на всех €чейках сетки
	//рассчитываетс€ по решению (UU) и наклонам (VV) на всех €чейках сетки
	//найденна€ скорость изменени€ решени€ и наклонов умножаетс€ на cft
	//(дл€ реализации €вного метода Ёйлера следует положить cft = ptrprm->tau)
	void step(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& DSOL, \
		const double cft);
};



#endif
