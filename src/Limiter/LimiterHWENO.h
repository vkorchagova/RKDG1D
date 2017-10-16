#ifndef LIMITERHWENO_H_
#define LIMITERHWENO_H_

#include <vector>

#include "Indicator.h"
#include "Limiter.h"

#define HWENO LimiterHWENO

using namespace std;

class LimiterHWENO :
    public Limiter
{
private:
    double wg; //—тепень, в которую возвод€тс€ слагаемые при расчете лимитированных наклонов
    
	const double weps = 1e-6; //ћалое положительное число

	//–≈јЋ»«ј÷»я ¬»–“”јЋ№Ќќ… ‘”Ќ ÷»» 
	//—обственно, расчет значений новых наклонов всех компонент решени€
	//заполнение значений новых наклонов (Vcorr[cell]) на конкретной (cell) €чейке сетки
	//рассчитываетс€ по решению (UU) и наклонам (VV) на всех €чейках сетки
    void CalculateBound(const vector<vector<double>>& UU, const vector<vector<double>>& VV, const int cell);

	double centrdiff(double yL, double yR, double h)
	{
		return 0.5 * (yR - yL) / h;
	}

	double diff2(double yL, double y, double yR, double h)
	{
		return (yR - 2.0*y + yL) / (h*h);
	}

	double backbackdiff(double yLL, double yL, double y, double h)
	{
		return 0.5 * (3.0*y - 4.0*yL + yLL) / h;
	}

	double forwforwdiff(double y, double yR, double yRR, double h)
	{
		return 0.5 * (-3.0*y + 4.0*yR - yRR) / h;
	}


public:
	// онструктор (prm - объект, содержащий параметры задачи, 
	//dimension - число консервативных переменных,
	//ind - объект, задающий индикаторную функцию,
	//degree - параметр лимитера, присваемый переменной wg)
    //prb - указатель на задачу,
    //LimiterHWENO(const BaseParams& prm, const int dimension, const Indicator& ind, double degree);
    LimiterHWENO(const BaseParams& prm, const Problem& prb, \
                 const Indicator& ind, double degree);

	//ƒеструктор
    ~LimiterHWENO();
};

#endif

