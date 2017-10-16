#ifndef LIMITERFINDIFF_H_
#define LIMITERFINDIFF_H_


#include "Limiter.h"

#define FinDiff LimiterFinDiff

class LimiterFinDiff :
    public Limiter
{
private:
	//–≈јЋ»«ј÷»я ¬»–“”јЋ№Ќќ… ‘”Ќ ÷»» 
	//—обственно, расчет значений новых наклонов всех компонент решени€
	//заполнение значений новых наклонов (Vcorr[cell]) на конкретной (cell) €чейке сетки
	//рассчитываетс€ по решению (UU) и наклонам (VV) на всех €чейках сетки
	//¬≈«ƒ≈ ”—“јЌј¬Ћ»¬јё“—я Ќ”Ћ≈¬џ≈ Ќј ЋќЌџ
	void CalculateBound(const vector<vector<vector<double>>>& SOL, \
		const int cell);
public:
	// онструктор (prm - объект, содержащий параметры задачи, 
	//dimension - число консервативных переменных,
	//ind - объект, задающий индикаторную функцию)
	//ƒолжен вызывать индикатор IndEverywhere, срабатывающий всюду!
    LimiterFinDiff(const BaseParams& prm, const Problem& prb, const Indicator& ind);

	//ƒеструктор
    ~LimiterFinDiff();
};

#endif
