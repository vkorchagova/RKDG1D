#ifndef FLUXHLL_H_
#define FLUXHLL_H_

#include <vector>
#include <cmath>
#include <algorithm>

#include "Flux.h"
#include "ProblemGas1D.h"

using namespace std;

#define HLL FluxHLL

class FluxHLL :
    public Flux
{
private:

    //—обственные числа на всей сетке (между i-й и (i-1)-й €чейками) 
    vector<vector<double>> L;

    //ѕотоки HLL на левом и правом конце €чейки
    vector<double> hllL, hllR;

    //¬ременные переменные: векторы - скачок решени€ на границах i-й €чейки
    vector<double> LdU, RdU;

	//—пособ вычислени€ скорости на границе €чеек
	SoundVelType SoundVel;

public:

    // онструктор (prm - объект, содержащий параметры задачи)
	FluxHLL(const BaseParams& prm, Problem& prb, SoundVelType soundvel);

    //ƒеструктор
    ~FluxHLL();

	void step(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& DSOL, \
		const double cft);

};


#endif