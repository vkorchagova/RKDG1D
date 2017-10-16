#ifndef INDICATORHARTEN_H_
#define INDICATORHARTEN_H_

#include "Indicator.h"
#include "defs.h"

#define Harten IndicatorHarten

class IndicatorHarten :
	public Indicator
{

private:
	double mult; //множитель дл€ сравнени€ наклонов в соседних €чейках

	//"„увствительна€" переменна€, по которой рассчитываетс€ значение индикаторной функции
	var sens;

public:
	// онструктор 
	//  prm - объект, содержащий параметры задачи, 
	//  val - "чувствительна€" переменна€
	//  alpha - коэффициент в индикаторе ’артена
	IndicatorHarten(const BaseParams& prm, const Problem& prb, const var val, const double alpha);

	//ƒеструктор
	~IndicatorHarten();

	//–≈јЋ»«ј÷»я ¬»–“”јЋ№Ќќ… ‘”Ќ ÷»» 
	//—обственно, расчет значений индикаторной функции
	//заполнение значений индикаторной функции (Ind) на всех €чейках сетки
	//рассчитываетс€ по решению (UU) и наклонам (VV) на всех €чейках сетки
	void calc_indicator(const vector<vector<double>>& UU, \
		const vector<vector<double>>& VV, \
		vector<double>& Ind) const;
};

#endif

