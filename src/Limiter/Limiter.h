#ifndef LIMITER_H_
#define LIMITER_H_

#include <vector>

#include "Indicator.h"
#include "Problem.h"

using namespace std;

class Limiter
{   
    private:
    protected:

		//¬ектор значений индикаторной функции во всех €чейках сетки
        vector<double> Ind;

		//—корректированные наклоны всех компонент решени€ во всех €чейках
		vector<vector<vector<double>>> SOLcorr;

		//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я 
		//—обственно, расчет значений новых наклонов всех компонент решени€
		//заполнение значений новых наклонов (Vcorr[cell]) на конкретной (cell) €чейке сетки
		//рассчитываетс€ по решению (UU) и наклонам (VV) на всех €чейках сетки
        virtual void CalculateBound(const vector<vector<vector<double>>>& SOL, \
            const int cell) = 0;
    
		//”казатель на объект, содержащий параметры задачи
		const BaseParams* ptrprm;

		//”казатель на объект, содержащий параметры задачи
		const Problem* ptrprb;

		//”казатель на индикатор (вызываетс€ дл€ расчета значений индикаторной функции)
		const Indicator* ptrInd;

    public:

		// онструктор (prm - объект, содержащий параметры задачи, 
		//dimension - число консервативных переменных,
		//ind - объект, задающий индикаторную функцию)
        Limiter(const BaseParams& prm, const Problem& prb, const Indicator& ind);
        
		//ƒеструктор
		~Limiter();

		//‘ункци€, вызывающа€ (при необходимости) расчет лимитированных значений
		//(заполнение массива VV)
		//рассчитываютс€ по решению (UU) и наклонам (VV) на всей сетке
		void Bound(vector<vector<vector<double>>>& SOL);
};

#endif
