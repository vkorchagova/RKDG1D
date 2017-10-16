#ifndef INDICATOR_H_
#define INDICATOR_H_

#include <vector>

#include "defs.h"
#include "Problem.h"

using namespace std;

class Indicator
{
private:

protected:
	//”казатель на объект с параметрами задачи
    const BaseParams* ptrprm;
	const Problem* ptrprb;
public:
	// онструктор (prm - объект, содержащий параметры задачи, prb - ссылка на решаемую задачу)
	Indicator(const BaseParams& prm, const Problem& prb);
    
	//ƒеструктор
	~Indicator();
    	
	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я 
	//—обственно, расчет значений индикаторной функции
	//заполнение значений индикаторной функции (Ind) на всех €чейках сетки
	//рассчитываетс€ по решению (UU) и наклонам (VV) на всех €чейках сетки
    virtual void calc_indicator(const vector<vector<vector<double>>>& SOL, \
        vector<double>& Ind) const = 0;
};

#endif
