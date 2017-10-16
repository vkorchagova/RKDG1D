#ifndef FLUXHLLE_H_
#define FLUXHLLE_H_

#include <vector>
#include <cmath>
#include <algorithm>

#include "Flux.h"
#include "ProblemGas1D.h"

using namespace std;

#define HLLE FluxHLLE

class FluxHLLE :
public Flux
{
private:
    
	//Собственные числа на всей сетке (между i-й и (i-1)-й ячейками) 
    vector<vector<double>> L;
    
	//Потоки HLL на левом и правом конце ячейки
    vector<double> hllL, hllR;
    
	//Временные переменные: векторы - скачок решения на границах i-й ячейки
    vector<double> LdU, RdU;
    
public:
    
	//Конструктор (prm - объект, содержащий параметры задачи)
    FluxHLLE(const BaseParams& prm, Problem& prb);
    
	//Деструктор
    ~FluxHLLE();
    
    void step(const vector<vector<double>>& UU, const vector<vector<double>>& VV, \
              vector<vector<double>>& DU, vector<vector<double>>& DV, \
              const double cft);
    
};


#endif
