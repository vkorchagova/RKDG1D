#ifndef LIMITERWENOH_H_
#define LIMITERWENOH_H_

#include <vector>

#include "Indicator.h"
#include "Limiter.h"

#define WENOH LimiterWENOH

using namespace std;

class LimiterWENOH :
    public Limiter
{
private:
    double wg; //—тепень, в которую возвод€тс€ слагаемые при расчете лимитированных наклонов
    
	const double weps = 1e-8; //ћалое положительное число

	//–≈јЋ»«ј÷»я ¬»–“”јЋ№Ќќ… ‘”Ќ ÷»» 
	//—обственно, расчет значений новых наклонов всех компонент решени€
	//заполнение значений новых наклонов (Vcorr[cell]) на конкретной (cell) €чейке сетки
	//рассчитываетс€ по решению (UU) и наклонам (VV) на всех €чейках сетки
    void CalculateBound(const vector<vector<double>>& UU, const vector<vector<double>>& VV, const int cell);
public:
	// онструктор (prm - объект, содержащий параметры задачи, 
	//dimension - число консервативных переменных,
	//ind - объект, задающий индикаторную функцию,
	//degree - параметр лимитера, присваемый переменной wg)
    //LimiterWENOH(const BaseParams& prm, const int dimension, const Indicator& ind, double degree);
    LimiterWENOH(const BaseParams& prm, const Problem& prb, \
                 const Indicator& ind, double degree);
	//ƒеструктор
    ~LimiterWENOH();
};

#endif

