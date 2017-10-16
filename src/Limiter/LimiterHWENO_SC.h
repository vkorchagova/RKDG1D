#ifndef LIMITERHWENO_SC_H_
#define LIMITERHWENO_SC_H_

#include <vector>

#include "Indicator.h"
#include "Limiter.h"

#define HWENO_SC LimiterHWENO_SC

using namespace std;

class LimiterHWENO_SC :
	public Limiter
{
private:
	double wg; //Степень, в которую возводятся слагаемые при расчете лимитированных моментов

	const double weps = 1e-6; //Малое положительное число

	//РЕАЛИЗАЦИЯ ВИРТУАЛЬНОЙ ФУНКЦИИ 
	//Собственно, расчет значений новых моментов всех компонент решения
	void CalculateBound(const vector<vector<vector<double>>>& SOL, const int cell);

public:
	//Конструктор (prm - объект, содержащий параметры задачи, 
	//dimension - число консервативных переменных,
	//ind - объект, задающий индикаторную функцию,
	//degree - параметр лимитера, присваемый переменной wg)
    LimiterHWENO_SC(const BaseParams& prm, const Problem& prb, const Indicator& ind, double degree);

	//Деструктор
	~LimiterHWENO_SC();
};

#endif
