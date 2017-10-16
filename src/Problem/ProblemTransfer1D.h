#ifndef PROBLEMTRANSFER1D_H_
#define PROBLEMTRANSFER1D_H_

#define Transfer1D ProblemTransfer1D

#include "Problem.h"
class ProblemTransfer1D :
	public Problem
{
private:
	
protected:

	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я
	//¬ычисление конвективного потока по заданному вектору решени€
	virtual void getFlux(const vector<double>& U, vector<double>& Flux) const;

public:
	ProblemTransfer1D(const BaseParams& prm, int order, vector<std::function<double(const double)>> initv, \
		std::function<double(const double)> vel);
	~ProblemTransfer1D();

	//‘ункци€, определ€юща€ скорость переноса
	std::function<double(const double)> physflux;

	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я
	//ѕредельное значение (sd = left/right) q-й компоненты решени€,
	//рассчитываемое по решению (U) и наклону (V)
	virtual inline double side_val(const vector<vector<double>>& sol, var q, side sd) const;
	
	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я
	//¬озвращает значение переменной q в точке x, вычисл€емое по заданным в объекте Param начальным услови€м
	// cft = 0 --- среднее значение
	// cft = 1 --- наклон
	virtual double initial_var(const double x, const var q, const int cft) const;

	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я
	//ѕечать в телефайл
	//virtual void more_to_file(ostream& str, \
	//	const vector<vector<double>>& UU, const vector<vector<double>>& VV, int cell) const;

	void EigenMatricies(const vector<vector<double>>& sol, \
		vector<vector<double>>& LL, \
		vector<vector<double>>& RR) const {};
};

#endif