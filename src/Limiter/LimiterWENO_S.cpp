#include "LimiterWENO_S.h"
#include "Integrator.h"


//КОНСТРУКТОР Лимитера WENO_S
LimiterWENO_S::LimiterWENO_S(const BaseParams& prm, const Problem& prb, \
    const Indicator& ind, double degree) : Limiter(prm, prb, ind)
{
	wg = degree;
}

//ДЕСТРУКТОР Лимитера
LimiterWENO_S::~LimiterWENO_S()
{}


void LimiterWENO_S::CalculateBound(const vector<vector<vector<double>>>& SOL, const int cell)
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
	double ak = 0.001;
	double gamma0 = ak;
	double gamma1 = 1.0 - 2.0*ak;
	double gamma2 = ak;

	if (nshape==3) 
	// Цикл по консервативным переменным 
	for (int val = 0; val < dim; ++val)
	{
		// Переобозначение моментов решения с ячеек шаблона.
		//const double& ul = leftSol[0][val];
		//const double& u = troublSol[0][val];
		//const double& ur = rightSol[0][val];
		const double& vl = leftSol[1][val];
		const double& v  = troublSol[1][val];
		const double& vr = rightSol[1][val];
		const double& wl = leftSol[2][val];
		const double& w = troublSol[2][val];
		const double& wr = rightSol[2][val];

		// Вычисляем индикаторы гладкости полиномов решения в проблемной ячейке. 
		double beta0 = 4.0 * ((vl + 6.0 * wl)*(vl + 6.0 * wl) + 39.0 * wl*wl);
		double beta1 = 4.0 * (v*v + 39.0 * w*w);
		double beta2 = 4.0 * ((vr - 6.0 * wr)*(vr - 6.0 * wr) + 39.0 * wr*wr);

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

		// Честно сделанный пересчёт моментов с учётом привязки базисных функций к ячейкам. 
		double t1, t2, t3,t4;
		
		t1 = w0*leftSol[2][val];
		t2 = w1*troublSol[2][val];
		t3 = w2*rightSol[2][val];
		t4 = t1 + t2 + t3;
		if (fabs(t4) < weps) t4 = 0;
		SOLcorr[cell][2][val] = t4;//(w0*leftSol[2][val] + w1*troublSol[2][val] + w2*rightSol[2][val]);
		
		t1 = w0*(leftSol[1][val] + 6.0 * leftSol[2][val]);
		t2 = w1*troublSol[1][val];
		t3 = w2*(rightSol[1][val] - 6.0 * rightSol[2][val]);
		t4 = t1 + t2 + t3;
		if (fabs(t4) < weps) t4 = 0;
		SOLcorr[cell][1][val] = t4;// (w0*(leftSol[1][val] + 6 * leftSol[2][val]) + w1*troublSol[1][val] + w2*(rightSol[1][val] - 6 * rightSol[2][val]));
			
	}
	else
		for (int val = 0; val < dim; ++val)
		{
			// Denoting slopes of the solution at the 3-point stencil.
			const double& vl = leftSol[1][val];
			const double& v = troublSol[1][val];
			const double& vr = rightSol[1][val];

			// Calculating the smothness indicators for solution at the troubled cell. 
			/*double beta0 = vl*vl *2.0 / h;
			double beta1 = v*v   *2.0 / h;
			double beta2 = vr*vr *2.0 / h;*/
			/*???*/
			double beta0 = 4.0*vl*vl;
			double beta1 = 4.0*v*v;
			double beta2 = 4.0*vr*vr;
			

			// Calculating the nonlinear weights. Normalization of weights will be carried out (applied?) during the final weightened interpolation.
			double w0 = gamma0 / pow(weps + beta0, wg);
			double w1 = gamma1 / pow(weps + beta1, wg);
			double w2 = gamma2 / pow(weps + beta2, wg);

			// The resulting slope of solution is evaluated as the linear combination of the existed slopes.
			SOLcorr[cell][1][val] = (w0*vl + w1*v + w2*vr) / (w0 + w1 + w2);
		}
}
