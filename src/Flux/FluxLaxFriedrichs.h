#ifndef FLUXLAXFRIEDRICHS_H_
#define FLUXLAXFRIEDRICHS_H_

#include <vector>
#include <cmath>
#include <algorithm>

#include "Flux.h"
#include "ProblemGas1D.h"
#include "ProblemMHD1D.h"

#define LaxFriedrichs FluxLaxFriedrichs

using namespace std;

class FluxLaxFriedrichs :
    public Flux
{
private:

    //—обственные числа на всей сетке (между i-й и (i-1)-й €чейками) 
    vector<vector<double>> L;

    //ѕотоки Ћакса-‘ридрихса на левом и правом конце €чейки
    vector<double> lfL, lfR;

    //¬ременные переменные: векторы - скачок решени€ на границах i-й €чейки
    vector<double> LdU, RdU;

	//—пособ вычислени€ скорости на границе €чеек
	SoundVelType SoundVel;

public:

    // онструктор  (prm - объект, содержащий параметры задачи)
	FluxLaxFriedrichs(const BaseParams& prm, Problem& prb, SoundVelType soundvel);

    //ƒеструктор
    ~FluxLaxFriedrichs();

	void step(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& DSOL, \
		const double cft);
};

#endif