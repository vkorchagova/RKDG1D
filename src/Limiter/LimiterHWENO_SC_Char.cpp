#include "LimiterHWENO_SC_Char.h"
#include "Integrator.h"


//КОНСТРУКТОР Лимитера HWENO_SC_Char
LimiterHWENO_SC_Char::LimiterHWENO_SC_Char(const BaseParams& prm, const Problem& prb, \
    const Indicator& ind, double degree) : Limiter(prm, prb, ind)
{
	wg = degree;
}

//ДЕСТРУКТОР Лимитера
LimiterHWENO_SC_Char::~LimiterHWENO_SC_Char()
{}


void LimiterHWENO_SC_Char::CalculateBound(const vector<vector<vector<double>>>& SOL, const int cell)
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

	//Матрицы для собственных векторов.
	vector<vector<double>> LL; 
	vector<vector<double>> RR;
	LL.resize(dim);
	RR.resize(dim);
	for (int j = 0; j < dim; ++j)
	{
		LL[j].resize(dim);
		RR[j].resize(dim);
	}// for j

	// Заполняем матрицы собственных векторов.
	ptrprb->EigenMatricies(troublSol, LL, RR);

	// Линейные веса.
	double ak = 0.001;
	double gamma0 = ak;
	double gamma1 = 1.0 - 2.0*ak;
	double gamma2 = ak;

	
	// Массив не- и лимитированных полиномов в характеристическом пространстве.
	vector<vector<vector<double>>> auxSOL;
	vector<vector<double>> auxSOLcorr;
	auxSOL.resize(3);
	auxSOLcorr.resize(nshape);
	for (int j = 0; j < nshape; ++j)
	{
		auxSOLcorr[j].resize(dim);
	}// for j
	for (int i = 0; i < 3; ++i)
	{
		auxSOL[i].resize(nshape);
		for (int j = 0; j < nshape; ++j)
		{
			auxSOL[i][j].resize(dim);
		}// for j		
	}// for i

	// Подготовка к расчёту по методу HWENO_SC --- проекция решения на базис собственных векторов.
	for (int j = 0; j < nshape; ++j)
	{
		prodMatrVec(LL, leftSol[j],   auxSOL[0][j]);
		prodMatrVec(LL, troublSol[j], auxSOL[1][j]);
		prodMatrVec(LL, rightSol[j],  auxSOL[2][j]);
	}


	// Цикл по консервативным переменным 
	for (int val = 0; val < dim; ++val)
	{
		// Переобозначение моментов решения с ячеек шаблона.
		const double& ul = auxSOL[0][0][val];
		const double& u = auxSOL[1][0][val];
		const double& ur = auxSOL[2][0][val];
		const double& vl = auxSOL[0][1][val];
		const double& v = auxSOL[1][1][val];
		const double& vr = auxSOL[2][1][val];
		
		// Создаём переменные, куда потом запишем всё необходимое. Формат вызван вилкой через if(nshape).
		double beta0, beta1, beta2;
		// Массив для боковых полиномов: интерполянов или просто модифицированных.
		vector<vector<double>> SIDE_auxSOL;
		SIDE_auxSOL.resize(nshape);
		for (int j = 0; j < nshape; ++j)
			SIDE_auxSOL[j].resize(2);

		double p;
		
		// Дальше считаем по методике в зависимости от nshape.
		if (nshape == 2)
		{
			// Evaluating modified with LSM polynomials. Operating with slopes seems to be enough.
			SIDE_auxSOL[1][0] = (2.0*vl / h + 12.0 * (u - ul) / h) / 13.0;
			p = 2.0 *v / h;
			SIDE_auxSOL[1][1] = (2.0*vr / h + 12.0 * (ur - u) / h) / 13.0;

			// Calculating the smothness indicators for modified polynomials at the troubled cell. 
			beta0 = SIDE_auxSOL[1][0]*SIDE_auxSOL[1][0] * h * h;
			beta1 = p*p * h * h;
			beta2 = SIDE_auxSOL[1][1]*SIDE_auxSOL[1][1] * h * h;

		}//if nshape
		
		if (nshape == 3)
		{
			const double& wl = auxSOL[0][2][val];
			const double& w = auxSOL[1][2][val];
			const double& wr = auxSOL[2][2][val];


			// Модифицированные полиномы соседних ячеек (по МНК)
			SIDE_auxSOL[0][0] = (u + 192.0*ul - 2.0*(vl + 3.0*wl)) / 193.0;
			SIDE_auxSOL[1][0] = (6.0*u - 6.0*ul + 181.0*vl - 36.0*wl) / 193.0;
			SIDE_auxSOL[2][0] = (30.0*u - 30.0*ul - 60.0*vl + 13.0*wl) / 193.0;

			SIDE_auxSOL[0][1] = (u + 2.0 * (96.0 * ur + vr - 3.0 * wr)) / 193.0;
			SIDE_auxSOL[1][1] = (-6.0*u + 6.0*ur + 181.0*vr + 36.0*wr) / 193.0;
			SIDE_auxSOL[2][1] = (30.0*u - 30.0*ur + 60.0*vr + 13.0*wr) / 193.0;


			// Вычисляем индикаторы гладкости полиномов решения в проблемной ячейке. 
			beta0 = 4.0 * ((SIDE_auxSOL[1][0] + 6.0 * SIDE_auxSOL[2][0])*(SIDE_auxSOL[1][0] + 6.0 * SIDE_auxSOL[2][0]) + 39.0 * SIDE_auxSOL[2][0]*SIDE_auxSOL[2][0]);
			beta1 = 4.0 * (v*v + 39.0 * w*w);
			beta2 = 4.0 * ((SIDE_auxSOL[1][1] - 6.0 * SIDE_auxSOL[2][1])*(SIDE_auxSOL[1][1] - 6.0 * SIDE_auxSOL[2][1]) + 39.0 * SIDE_auxSOL[2][1]*SIDE_auxSOL[2][1]);
		}// if nshape
		// Вычисляем нелинейные веса. Нормировка будет осуществлена ниже.
		double w0 = gamma0 / pow(weps + beta0, wg);
		double w1 = gamma1 / pow(weps + beta1, wg);
		double w2 = gamma2 / pow(weps + beta2, wg);
		// Нормировка весов.
		double WWW = (w0 + w1 + w2);
		w0 = w0 / WWW; w1 = w1 / WWW; w2 = w2 / WWW;
		
		
		// Честно сделанный пересчёт моментов с учётом привязки базисных функций к ячейкам. 

		double t1, t2, t3,t4;
		t1 = w0*SIDE_auxSOL[nshape-1][0];
		t2 = w1*troublSol[nshape-1][val];
		t3 = w2*SIDE_auxSOL[nshape-1][1];
		t4 = t1 + t2 + t3;
		//if (fabs(t4) < weps) t4 = 0;
		auxSOLcorr[nshape-1][val] = t4;//(w0*leftSol[2][val] + w1*troublSol[2][val] + w2*rightSol[2][val]);
		
		if (nshape == 3)
		{
			t1 = w0*(SIDE_auxSOL[1][0] + 6.0 * SIDE_auxSOL[2][0]);
			t2 = w1*troublSol[1][val];
			t3 = w2*(SIDE_auxSOL[1][1] - 6.0 * SIDE_auxSOL[2][1]);
			t4 = t1 + t2 + t3;
			//if (fabs(t4) < weps) t4 = 0;
			auxSOLcorr[1][val] = t4;// (w0*(leftSol[1][val] + 6 * leftSol[2][val]) + w1*troublSol[1][val] + w2*(rightSol[1][val] - 6 * rightSol[2][val]));
		}//if nshapae

	}// for val

	// Проектируем лимитированный полином обратно в пространство консервативных переменных.
	for (int j = 1; j < nshape; ++j)
		prodMatrVec(RR, auxSOLcorr[j], SOLcorr[cell][j]);
}
