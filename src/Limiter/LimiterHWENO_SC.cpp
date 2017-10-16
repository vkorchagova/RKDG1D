#include "LimiterHWENO_SC.h"
#include "Integrator.h"


//КОНСТРУКТОР Лимитера HWENO_SC
LimiterHWENO_SC::LimiterHWENO_SC(const BaseParams& prm, const Problem& prb, \
    const Indicator& ind, double degree) : Limiter(prm, prb, ind)
{
	wg = degree;
}

//ДЕСТРУКТОР Лимитера
LimiterHWENO_SC::~LimiterHWENO_SC()
{}


void LimiterHWENO_SC::CalculateBound(const vector<vector<vector<double>>>& SOL, const int cell)
{
	// Некоторые базовые параметры.
	int nx = ptrprm->nx;
	double h = ptrprm->h;
	int dim = ptrprb->dim;
	int nshape = ptrprb->nshape;

	//Определяем соседей проблемной ячейки. При необходимости, посредством фиктивных ячеек.
	const vector<vector<double>>& leftSol = (cell > 0) ? SOL[cell - 1] : SOL[nx];
	const vector<vector<double>>& troublSol = SOL[cell];
	const vector<vector<double>>& rightSol = (cell < nx - 1) ? SOL[cell + 1] : SOL[nx + 1];

	// Линейные веса.
	double ak = 0.01;
	double gamma0 = ak;
	double gamma1 = 1.0 - 2.0*ak;
	double gamma2 = ak;

	// Цикл по консервативным переменным 
	for (int val = 0; val < dim; ++val)
	{
		// Переобозначение моментов решения с ячеек шаблона.
		const double& ul = leftSol[0][val];
		const double& u = troublSol[0][val];
		const double& ur = rightSol[0][val];
		const double& vl = leftSol[1][val];
		const double& v  = troublSol[1][val];
		const double& vr = rightSol[1][val];
		const double& wl = leftSol[2][val];
		const double& w = troublSol[2][val];
		const double& wr = rightSol[2][val];

		// Модифицированные полиномы соседних ячеек (по МНК)
		double pLu = (u + 192.0*ul - 2.0*(vl + 3.0*wl)) / 193.0;
		double pLv = (6.0*u - 6.0*ul + 181.0*vl - 36.0*wl) / 193.0;
		double pLw = (30.0*u - 30.0*ul - 60.0*vl + 13.0*wl) / 193.0;
		
		double pRu = (u + 2 * (96 * ur + vr - 3 * wr)) / 193.0;
		double pRv = (-6.0*u + 6.0*ur + 181.0*vr + 36.0*wr) / 193.0;
		double pRw = (30.0*u - 30.0*ur + 60.0*vr + 13.0*wr) / 193.0;


		// Вычисляем индикаторы гладкости полиномов решения в проблемной ячейке. 
		double beta0 = 4.0 * ((pLv + 6.0 * pLw)*(pLv + 6.0 * pLw) + 39.0 * pLw*pLw);
		double beta1 = 4.0 * (v*v + 39.0 * w*w);
		double beta2 = 4.0 * ((pRv - 6.0 * pRw)*(pRv - 6.0 * pRw) + 39.0 * pRw*pRw);

		// Вычисляем нелинейные веса. Нормировка будет осуществлена ниже.
		double w0 = gamma0 / pow(weps + beta0, wg);
		double w1 = gamma1 / pow(weps + beta1, wg);
		double w2 = gamma2 / pow(weps + beta2, wg);
		// Нормировка весов.
		double WWW = (w0 + w1 + w2);
		w0 = w0 / WWW; w1 = w1 / WWW; w2 = w2 / WWW;
		
		//// Эксперимент против логики: пробуем игнорировать привязку базисных функций к ячейкам. 
		//for (int shape = 1; shape < ptrprb->nshape; shape++)
		//	SOLcorr[cell][shape][val] = (w0*leftSol[shape][val] + w1*troublSol[shape][val] + w2*rightSol[shape][val]) / (w0 + w1 + w2);

		//// Честно сделанный пересчёт моментов с учётом привязки базисных функций к ячейкам. 
		double t1, t2, t3, t4;
		
		if (ptrprb->nshape > 2)
		{
			t1 = w0*pLw;
			t2 = w1*troublSol[2][val];
			t3 = w2*pRw;
			t4 = t1 + t2 + t3;
			//if (fabs(t4) < weps) t4 = 0;
			SOLcorr[cell][2][val] = t4;//(w0*leftSol[2][val] + w1*troublSol[2][val] + w2*rightSol[2][val]);
		}

		t1 = w0*(pLv + 6.0 * pLw);
		t2 = w1*troublSol[1][val];
		t3 = w2*(pRv - 6.0 * pRw);
		t4 = t1 + t2 + t3;
		//if (fabs(t4) < weps) t4 = 0;
		SOLcorr[cell][1][val] = t4;// (w0*(leftSol[1][val] + 6 * leftSol[2][val]) + w1*troublSol[1][val] + w2*(rightSol[1][val] - 6 * rightSol[2][val]));
			
	}
}
