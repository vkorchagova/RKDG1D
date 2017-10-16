#include "LimiterWENOH.h"



//КОНСТРУКТОР Лимитера WENOH
LimiterWENOH::LimiterWENOH(const BaseParams& prm, const Problem& prb,\
    const Indicator& ind, double degree) : Limiter(prm, prb, ind)
{
    wg = degree;
}
//ДЕСТРУКТОР Лимитера
LimiterWENOH::~LimiterWENOH()
{};


void LimiterWENOH::CalculateBound(const vector<vector<double>>& UU, const vector<vector<double>>& VV, const int cell)
{
    double m0, m1, m2, m3, m4, znam;

    vector<double> pows(5); //это - от WENOH, с размерностью задачи не связано!
   
	int nx = ptrprm->nx;
	double h = ptrprm->h;

	const vector<double>& leftV  = (cell > 0) ? VV[cell - 1] : VV[nx];
    const vector<double>& rightV = (cell < nx - 1) ? VV[cell + 1] : VV[nx + 1];

	const vector<double>& leftU = (cell > 0) ? UU[cell - 1] : UU[nx];
	const vector<double>& rightU = (cell < nx - 1) ? UU[cell + 1] : UU[nx + 1];

    for (int val = 0; val < dim; ++val)
    {
        m0 = 2.0 / h * VV[cell][val];
        m1 = 2.0 / h * leftV[val];
        m2 = 2.0 / h * rightV[val];
		m3 = (UU[cell][val] - leftU[val]) / h;
		m4 = (rightU[val] - UU[cell][val]) / h;

        pows[0] = pow((weps + fabs(m0)), -wg);
        pows[1] = pow((weps + fabs(m1)), -wg);
        pows[2] = pow((weps + fabs(m2)), -wg);
        pows[3] = pow((weps + fabs(m3)), -wg);
        pows[4] = pow((weps + fabs(m4)), -wg);

        znam = pows[0] + pows[1] + pows[2] + pows[3] + pows[4];

        Vcorr[cell][val] = 0.5 * h * (pows[0] * m0 + pows[1] * m1 + pows[2] * m2 + pows[3] * m3 + pows[4] * m4) / znam;
    }
}
