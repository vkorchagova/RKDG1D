#ifndef PROBLEM_H_
#define PROBLEM_H_

#include "defs.h"

class Problem
{
protected:
	//”казатель на объект с параметрами задачи
	const BaseParams* ptrprm;

public:
	//–азмерность задачи (количество консервативных переменных)
	int dim;

	// оличество базисных функций
	int nshape;

	function<double(double)> shapefunc[3];
	double shapefuncnorm2[3];

	//”казатели на функции, содержащие начальные услови€ (U) либо (\rho, Ux, Uy, Uz, p)
	vector<std::function<double(const double)>> init; 

	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я
	//ѕредельное значение (sd = left/right) q-й компоненты решени€,
	//рассчитываемое по решению (U) и наклону (V)
	virtual inline double side_val(const vector<vector<double>>& sol, var q, side sd) const = 0;
	virtual inline vector<double> gauss_val(const vector<vector<double>>& sol, double Lcoord) const = 0;


	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я
	//¬ычисление аналитического потока Flux по заданному решению U
	virtual void getFlux(const vector<double>& U, vector<double>& Flux) const = 0;
	virtual vector<double> getFlux(const vector<double>& U) const = 0;

	//јналитические потоки, вычисленные по решению:
	//  - в центре €чейки (FLUX), 
	//  - на левом конце (LFLUX), 
	//  - на правом конце (RFLUX)
	// ¬ позици€х n и (n+1) - фиктивные €чейки
	vector<vector<double>> FLUX, LFLUX, RFLUX;

	//«аполнение вектора потоков дл€ всех €чеек:
	//  - по решению в центах €чеек: conv_flux(UU)
	//  - по решению в центах €чеек, на левом и правом кра€х: conv_flux(UU, VV) с учетом фиктивных €чеек
//	void convFlux(const vector<vector<double>>& UU);
	void convFlux(const vector<vector<vector<double>>>& SOL);
		

	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я
	//¬озвращает вектор конс. переменных в точке x, вычисл€емый по заданным в объекте Param начальным услови€м
	virtual vector<double> initial_var(const double x) const = 0;


	virtual void EigenMatricies(const vector<vector<double>>& sol, \
		vector<vector<double>>& LL, \
		vector<vector<double>>& RR) const = 0;

	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я — ѕ”—“ќ… –≈јЋ»«ј÷»≈… ѕќ ”ћќЋ„јЌ»ё
	//ѕечать в телефайл
	virtual void more_to_file(ostream& str, \
		const vector<vector<vector<double>>>& SOL, int cell) const {	};


	Problem(const BaseParams& prm, int dimension, int nshapefunctions, vector<std::function<double(const double)>> initv);
	~Problem();
};

#endif