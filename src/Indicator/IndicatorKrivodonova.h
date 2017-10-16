#ifndef INDIKATORKRIVODONOVA_H_
#define INDICATORKRIVODONOVA_H_

#include "Indicator.h"
#include "defs.h"


#define Krivodonova IndicatorKrivodonova

class IndicatorKrivodonova :
    public Indicator
{
private:
	//"„увствительна€" переменна€, по которой рассчитываетс€ значение индикаторной функции
	var sens;

public:
	// онструктор (prm - объект, содержащий параметры задачи, val - "чувствительна€" переменна€
	//             prb - решаема€ задача (т.к. важно направление переноса!)
	IndicatorKrivodonova(const BaseParams& prm, const Problem& prb, const var val);

	//ƒеструктор
    ~IndicatorKrivodonova();
    
	//–≈јЋ»«ј÷»я ¬»–“”јЋ№Ќќ… ‘”Ќ ÷»» 
	//—обственно, расчет значений индикаторной функции
	//заполнение значений индикаторной функции (Ind) на всех €чейках сетки
	//рассчитываетс€ по решению (UU) и наклонам (VV) на всех €чейках сетки
	void calc_indicator(const vector<vector<vector<double>>>& SOL, \
		vector<double>& Ind) const;
};

#endif