#ifndef FLUX_H_
#define FLUX_H_

#include <vector>

#include "defs.h"
#include "Problem.h"
#include "ProblemGas1D.h"
#include "ProblemMHD1D.h"

using namespace std;

class Flux
{
private:
	//”казатели на векторы конв. потоков на трех €чейках
	vector<double> *ptr_leftfluxR, *ptr_myflux, *ptr_myfluxL, *ptr_myfluxR, *ptr_rightfluxL;

	//”казатели на векторы решени€ в трех €чейках
	const vector<vector<double>> *ptr_mysol, *ptr_leftsol, *ptr_rightsol;


protected:
	//”казатель на объект с общими параметрами
	const BaseParams* ptrprm;
			
	//”казатель на объект задача
	Problem* ptrprb;
	ProblemGas1D* ptrprb_toGas;
	ProblemMHD1D* ptrprb_toMHD;

	//ќбертка дл€ вызова одноименной функции из класса Problem
    double side_val(const vector<vector<double>>& sol, var q, side sd) const
	{
		return ptrprb->side_val(sol, q, sd);
	}

    vector<double> gauss_val(const vector<vector<double>>& sol, double Lcoord) const
	{
		return ptrprb->gauss_val(sol, Lcoord);
	}

	//—сылки на векторы конв. потоков на трех €чейках
	vector<double>& myfluxL() { return *ptr_myfluxL; };
	vector<double>& myfluxR() { return *ptr_myfluxR; };
	vector<double>& myflux()  { return *ptr_myflux; };
	vector<double>& leftfluxR()  { return *ptr_leftfluxR; };
	vector<double>& rightfluxL() { return *ptr_rightfluxL; };

	//—сылки на векторы решени€ на трех €чейках
	const vector<vector<double>>& mysol() { return *ptr_mysol; };
	const vector<vector<double>>& leftsol()  { return *ptr_leftsol; };
	const vector<vector<double>>& rightsol() { return *ptr_rightsol; };
	

	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я — –≈јЋ»«ј÷»≈…
	//”становка приватных указателей (чтобы работали вышеперечисленные ссылки): 
	// - на конв. потоки в своей €чейке
	// - на предельные (слева и справа) потоки в фиктивных €чейках
	// - на решение в своей €чейке и в соседних €чейках
	virtual void setlocsolflux(const vector<vector<vector<double>>>& SOL, const int cell);




public: 

	// онструктор (prm - объект, содержащий параметры задачи)
	Flux(const BaseParams& prm, Problem& prb);
	
	//ƒеструктор
	~Flux();

	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я 
	//«аполнение приращений решений (DU) и наклонов (DV) на всех €чейках сетки
	//рассчитываетс€ по решению (UU) и наклонам (VV) на всех €чейках сетки
	//найденна€ скорость изменени€ решени€ и наклонов умножаетс€ на cft
	//(дл€ реализации €вного метода Ёйлера следует положить cft = ptrprm->tau)
	virtual void step(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& DSOL, \
		const double cft) = 0;
};

#endif

