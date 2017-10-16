#include "LimiterHWENO.h"



//КОНСТРУКТОР Лимитера HWENO
//LimiterHWENO::LimiterHWENO(const BaseParams& prm, const int dimension, \
//    const Indicator& ind, double degree) : Limiter(prm, dimension, ind)
LimiterHWENO::LimiterHWENO(const BaseParams& prm, const Problem& prb,\
    const Indicator& ind, double degree) : Limiter(prm, prb, ind)
{
    wg = degree;
}

//ДЕСТРУКТОР Лимитера
LimiterHWENO::~LimiterHWENO()
{}


void LimiterHWENO::CalculateBound(const vector<vector<double>>& UU, const vector<vector<double>>& VV, const int cell)
{
	int nx = ptrprm->nx;
	double h = ptrprm->h;

	const vector<double>& leftU = (cell > 0) ? UU[cell - 1] : UU[nx];
	const vector<double>& rightU = (cell < nx - 1) ? UU[cell + 1] : UU[nx + 1];

    const vector<double>& leftV = (cell > 0) ? VV[cell - 1] : VV[nx];
    const vector<double>& rightV = (cell < nx - 1) ? VV[cell + 1] : VV[nx + 1];


// --------- TODO: correct linear weights
    double gamma0m = 0.09;//0.3;//
    double gamma1m = 0.9;//0.6; //
    double gamma2m = 0.01;//0.1; //

    double gamma0p = 0.01;//0.1;//
    double gamma1p = 0.9;//0.6; //
    double gamma2p = 0.09;//0.3;//

    int dim = ptrprb->dim;

	for (int val = 0; val < dim; ++val)
	{
		const double& yl  = leftU[val];
        const double& dyl = 2.0*leftV[val]/h;
		const double& y   = UU[cell][val];
        const double& dy   = 2.0*VV[cell][val]/h;
		const double& yr  = rightU[val];
        const double& dyr  = 2.0*rightV[val]/h;
		
        double p0m = (h*dyl + 3.0*yl + y)/4.0;
        double p0p = (-3.0*dyl*h - 5.0*yl + 9.0*y)/4.0;

		double p1m = 0.125*(3.0*yl + 6.0*y - yr);
		double p1p = 0.125*(-yl + 6.0*y + 3.0*yr);

        double p2m = (3.0*dyr*h-5.0*yr+9.0*y)/4.0;
        double p2p = (-dyr*h+3.0*yr+y)/4.0;
		
        double beta0 = (16.0*dyl*dyl*h*h+25.0*(y-yl)*(y-yl)+38.0*dyl*h*(yl-y))/3.0; //(4.0*yll*yll - 19.0*yll*yl + 25.0*yl*yl + 11.0*yll*y - 31.0*yl*y + 10.0*y*y) / 3.0;
		double beta1 = (4.0*yl*yl - 13.0*yl*y + 13.0*y*y + 5.0*yl*yr - 13.0*y*yr + 4.0*yr*yr) / 3.0;
        double beta2 = (16.0*dyr*dyr*h*h+25.0*(y-yr)*(y-yr)+38.0*dyr*h*(y-yr))/3.0;//(10.0*y*y - 31.0*y*yr + 11.0*y*yrr + 25.0*yr*yr - 19.0*yr*yrr + 4.0*yrr*yrr) / 3.0;
		
		double w0m = gamma0m / pow(weps + beta0, wg);
		double w1m = gamma1m / pow(weps + beta1, wg);
		double w2m = gamma2m / pow(weps + beta2, wg);
		
		double w0p = gamma0p / pow(weps + beta0, wg);
		double w1p = gamma1p / pow(weps + beta1, wg);
		double w2p = gamma2p / pow(weps + beta2, wg);

		double Uplus =  (p0p*w0p + p1p*w1p + p2p*w2p) / (w0p + w1p + w2p);
		double Uminus = (p0m*w0m + p1m*w1m + p2m*w2m) / (w0m + w1m + w2m);

        double ttt = 0.5 * (Uplus - Uminus);

        for (int shape = 1; shape < ptrprb->nshape; ++shape)
              SOLcorr[cell][shape][val] = ttt;
	}
}
