#ifndef TIMESTEP_H_
#define TIMESTEP_H_

#include <vector>

#include "Flux.h"
#include "Limiter.h"
#include "Boundary.h"

class Timestep
{
protected:
	//”казатель на объект с параметрами задачи
	const BaseParams* ptrprm;
	
	// оличество консервативных переменных
	int dim;

	//”казатель на √”
	const Boundary* ptrbnd;

	//¬екторы приращений средних значенй и наклонов:
	vector<vector<vector<double>>> DSOL;
	
	//—нос решени€ U, V на фиктивную €чейку	
	void ApplyBoundary(const vector<vector<vector<double>>>& SOL);

public:
	// онструктор
	Timestep(const BaseParams& prm, int dimension, int nshape, const Boundary& bnd);
	
	//ƒеструктор
	virtual ~Timestep();

	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я 
	//—обственно, шаг расчета соответствующим численным методом
	//¬ычисление нового решени€ (Unew) и наклонов (Vnew) на всех €чейках сетки
	//рассчитываетс€ по решению (U) и наклонам (V) на всех €чейках сетки.
	//Ўаг по времени tau
	//ƒл€ расчета потоков используетс€ численный поток method, лимитер lim.
	//(в реализации дл€ расчета приращений следует использовать method.step(...)
	virtual void runstep(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& SOLnew, \
		const double tau, Flux& method, Limiter& lim) = 0;	
};

#endif