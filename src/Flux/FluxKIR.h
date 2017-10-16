#ifndef FLUXKIR_H_
#define FLUXKIR_H_

#include <vector>
#include <cmath>
#include <algorithm>

#include "Flux.h"
#include "ProblemGas1D.h"

#define KIR FluxKIR

using namespace std;


class FluxKIR : 
	public Flux
{
private:
	//Ћевые и правые собственные векторы на всей сетке (между i-й и (i-1)-й €чейками)
    vector<vector<vector<double>>> LW, RW; 

	//—обственные числа на всей сетке (между i-й и (i-1)-й €чейками) 
    vector<vector<double>> L;    

	//¬ременные переменные: матрицы - результат перемножени€ RW * L * LW
    vector<vector<double>> LMatr, RMatr;

	//¬ременные переменные: векторы - скачок решени€ на границах i-й €чейки
	vector<double> LdU, RdU;

	//¬ременные переменные: векторы - результат перемножени€ Matr * dU
    vector<double> ProdL, ProdR;	

public:
	// онструктор  (prm - объект, содержащий параметры задачи)
	FluxKIR(const BaseParams& prm, Problem& prb);

	//ƒеструктор
    ~FluxKIR();

	//–≈јЋ»«ј÷»я ¬»–“”јЋ№Ќќ… ‘”Ќ ÷»» 
	//—обственно, шаг расчета
	//заполнение приращений решений (DU) и наклонов (DV) на всех €чейках сетки
	//рассчитываетс€ по решению (UU) и наклонам (VV) на всех €чейках сетки
	//найденна€ скорость изменени€ решени€ и наклонов умножаетс€ на cft
	//(дл€ реализации €вного метода Ёйлера следует положить cft = ptrprm->tau)
    void step(const vector<vector<double>>& UU, const vector<vector<double>>& VV, \
		vector<vector<double>>& DU, vector<vector<double>>& DV, \
        const double cft);
};

#endif
