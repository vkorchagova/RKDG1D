#ifndef FLUXHLLC_H_
#define FLUXHLLC_H_

#include <vector>
#include <cmath>
#include <algorithm>

#include "Flux.h"
#include "ProblemGas1D.h"

using namespace std;

#define HLLC FluxHLLC

class FluxHLLC :
    public Flux
{
private:

    //Собственные числа на всей сетке (между i-й и (i-1)-й ячейками) 
    vector<vector<double>> L;

	//Способ вычисления скорости на границе ячеек
	SoundVelType SoundVel;

    //Потоки HLLС на левом и правом конце ячейки
    vector<double> hllcL, hllcR;

    //HLLC-решения "со звездочкой"
    vector<double> UastL, UastR;

    //Временные переменные: векторы - HLLC-скачок решения на границах i-й ячейки
    vector<double> LdUast, RdUast;

    void getUast(const int cell, const double Sast, \
		const vector<vector<double>>& leftsol, \
		const vector<vector<double>>& mysol, \
        vector<double>& Uast);

    double getSast(const int cell, \
		const vector<vector<double>>& leftsol, \
		const vector<vector<double>>& mysol);
public:

    //Конструктор  (prm - объект, содержащий параметры задачи)
	FluxHLLC(const BaseParams& prm, Problem& prb, SoundVelType soundvel);

    //Деструктор
    ~FluxHLLC();

	void step(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& DSOL, \
		const double cft);

};

#endif