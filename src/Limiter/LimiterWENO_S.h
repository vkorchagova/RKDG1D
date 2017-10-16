#ifndef LIMITERWENO_S_H_
#define LIMITERWENO_S_H_

#include <vector>

#include "Indicator.h"
#include "Limiter.h"

#define WENO_S LimiterWENO_S

using namespace std;

class LimiterWENO_S :
	public Limiter
{
private:
	double wg; //Степень, в которую возводятся слагаемые при расчете лимитированных моментов

	const double weps = 1e-7; //Малое положительное число

	//РЕАЛИЗАЦИЯ ВИРТУАЛЬНОЙ ФУНКЦИИ 
	//Собственно, расчет значений новых моментов всех компонент решения
	void CalculateBound(const vector<vector<vector<double>>>& SOL, const int cell);

public:
	//Конструктор (prm - объект, содержащий параметры задачи, 
	//dimension - число консервативных переменных,
	//ind - объект, задающий индикаторную функцию,
	//degree - параметр лимитера, присваемый переменной wg)
    LimiterWENO_S(const BaseParams& prm, const Problem& prb, const Indicator& ind, double degree);

	//Деструктор
	~LimiterWENO_S();
};

#endif
