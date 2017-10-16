#include "ProblemMHD1D.h"

const int dimMHD1D = 7;

ProblemMHD1D::ProblemMHD1D(const BaseParams& prm, int order, double gam, double hx, \
	vector<std::function<double(const double)>> initv) : Problem(prm, dimMHD1D, order, initv)
{
	gamma = gam;
	Hx = hx;
}


ProblemMHD1D::~ProblemMHD1D()
{
}


void ProblemMHD1D::getFlux(const vector<double>& U, vector<double>& Flux) const
{
	double r = U[var::r];
	double rvx = U[var::rvx];
	double rvy = U[var::rvy];
	double rvz = U[var::rvz];
	double e = U[var::e];
	double Hy = U[var::Hy];
	double Hz = U[var::Hz];

	double vx = rvx / r;
	double vy = rvy / r;
	double vz = rvz / r;
	double vv = (vx*vx + vy*vy + vz*vz);
	double HH = (Hx*Hx + Hy*Hy + Hz*Hz);
	double press = (gamma - 1)*(e - 0.5*r*vv - 0.5*HH);
	double pressTot = press + 0.5*HH;



	Flux[var::r] = rvx;
	Flux[var::rvx] = rvx*vx - Hx*Hx + pressTot;
	Flux[var::rvy] = rvy*vx - Hx*Hy;
	Flux[var::rvz] = rvz*vx - Hx*Hz;
	Flux[var::e] = (e + pressTot)*vx - Hx*(Hx*vx + Hy*vy + Hz*vz);
	Flux[var::Hy] = Hy*vx - vy*Hx;
	Flux[var::Hz] = Hz*vx - vz*Hx;
	
}

vector<double> ProblemMHD1D::getFlux(const vector<double>& U) const
{
	vector<double> F(dim);
	
	double r = U[var::r];
	double rvx = U[var::rvx];
	double rvy = U[var::rvy];
	double rvz = U[var::rvz];
	double e = U[var::e];
	double Hy = U[var::Hy];
	double Hz = U[var::Hz];

	double vx = rvx / r;
	double vy = rvy / r;
	double vz = rvz / r;
	double vv = (vx*vx + vy*vy + vz*vz);
	double HH = (Hx*Hx + Hy*Hy + Hz*Hz);
	double press = (gamma - 1)*(e - 0.5*r*vv - 0.5*HH);
	double pressTot = press + 0.5*HH;



	F[var::r] = rvx;
	F[var::rvx] = rvx*vx - Hx*Hx + pressTot;
	F[var::rvy] = rvy*vx - Hx*Hy;
	F[var::rvz] = rvz*vx - Hx*Hz;
	F[var::e] = (e + pressTot)*vx - Hx*(Hx*vx + Hy*vy + Hz*vz);
	F[var::Hy] = Hy*vx - vy*Hx;
	F[var::Hz] = Hz*vx - vz*Hx;

	return F;
}


//«начени€ всех переменных на левой и правой границах €чеек
inline double ProblemMHD1D::side_val(const vector<vector<double>>& sol, var q, side sd) const
{
	double sgn = (sd == side::left) ? -1.0 : 1.0;
	double res = 0.0;

	switch (q)
	{
	case r:
	case rvx:
	case rvy:
	case rvz:
	case e:
	case Hy:
	case Hz:
		res = sol[0][q] + sgn * sol[1][q];
		if (nshape>2) 
			res += sol[2][q];
		break;

	case vx:
	case vy:
	case vz:
		res = (sol[0][q - 100] + sgn * sol[1][q - 100]) / side_val(sol, var::r, sd);
		if (nshape > 2)
			res += sol[2][q - 100] / side_val(sol, var::r, sd);
		break;

	case p:
		res = (gamma - 1) * (side_val(sol, var::e, sd) - \
			0.5*(sqr(side_val(sol, var::rvx, sd)) + sqr(side_val(sol, var::rvy, sd)) + \
				 sqr(side_val(sol, var::rvz, sd))) / side_val(sol, var::r, sd) - \
			0.5*(Hx*Hx + sqr(side_val(sol, var::Hy, sd) + sqr(side_val(sol, var::Hz, sd)))));

		break;
	}//switch
	return res;
}

//«начени€ всех переменных на €чейке
inline double ProblemMHD1D::val(const vector<vector<double>>& sol, var q) const
{
	double res = 0.0;

	switch (q)
	{
	case r:
	case rvx:
	case rvy:
	case rvz:
	case e:
	case Hy:
	case Hz:
		res = sol[0][q];
		break;

	case vx:
	case vy:
	case vz:
		res = (sol[0][q - 100]) / val(sol, var::r);
		break;

	case p:
		res = (gamma - 1) * (val(sol, var::e) - \
			0.5*(sqr(val(sol, var::rvx)) + sqr(val(sol, var::rvy)) + \
			sqr(val(sol, var::rvz))) / val(sol, var::r) - \
			0.5*(Hx*Hx + sqr(val(sol, var::Hy)) + sqr(val(sol, var::Hz))));
		break;
	}//switch
	return res;
}

inline vector<double> ProblemMHD1D::gauss_val(const vector<vector<double>>& sol, double Lcoord) const
{	
	vector<double> res(sol[0]); 

	if (nshape > 1)
		res += sol[1] * Lcoord;

	if (nshape > 2)
		res += sol[2] * (1.5*Lcoord*Lcoord - 0.5);

	return res;	
}


//скорость звука на границе (left)-й и (right)-й €чеек
double ProblemMHD1D::c_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const
{
	//return 0.5*(sqrt(ptrprm->gamma*side_val(UU[cell], VV[cell], var::p, side::left) / side_val(UU[cell], VV[cell], var::r, side::left)) + \
			//	sqrt(ptrprm->gamma*side_val(UU[cell - 1], VV[cell - 1], var::p, side::right) / side_val(UU[cell - 1], VV[cell - 1], var::r, side::right)));
	double h = h_av(solleft, solright);
	double v2 = v2_av(solleft, solright);
	return sqrt((gamma - 1.0)*(h - 0.5*v2));
	//return 0.5*(sqrt(ptrprm->gamma * pl / rhol) + sqrt(ptrprm->gamma * pr / rhor));
	//return sqrt(ptrprm->gamma * (pl + pr) / (rhol + rhor));        
}     
//скорость звука на €чейке
double ProblemMHD1D::c(const vector<vector<double>>& sol) const
{
	/*double hh = h(sol);
	double vv2 = v2(sol);														   „“ќ «ƒ≈—№ ѕ–ќ»—’ќƒ»“?!?!
	return sqrt((gamma - 1.0)*(hh - 0.5*vv2)); */     
	return sqrt( (val(sol, var::p) * gamma) / val(sol, var::r) );
}

//?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
//энтальни€ на границе (left)-й и (right)-й €чеек																				
double ProblemMHD1D::h_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const
{
	double rleft = side_val(solleft, var::r, side::right);
	double rright = side_val(solright, var::r, side::left);
	double hleft = (side_val(solleft, var::p, side::right) \
		+ side_val(solleft, var::e, side::right)) / rleft;
	double hright = (side_val(solright, var::p, side::left) \
		+ side_val(solright, var::e, side::left)) / rright;
	double sqrtrleft = sqrt(rleft);
	double sqrtrright = sqrt(rright);
	return (sqrtrleft*hleft + sqrtrright*hright) / (sqrtrleft + sqrtrright);
	//return 0.5*((side_val(Uright, Vright, var::p, side::left) + side_val(Uright, Vright, var::e, side::left)) / side_val(Uright, Vright, var::r, side::left) +
	//(side_val(Uleft, Vleft, var::p, side::right) + side_val(Uleft, Vleft, var::e, side::right)) / side_val(Uleft, Vleft, var::r, side::right));
}		
//энтальни€ на €чейке
double ProblemMHD1D::h(const vector<vector<double>>& sol) const
{
	return (val(sol, var::p) + val(sol, var::e)) / val(sol, var::r);
}

//скорость на границе €чеек
vector<double> ProblemMHD1D::v_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const
{
	double sqrtrleft = sqrt(side_val(solleft, var::r, side::right));
	double sqrtrright = sqrt(side_val(solright, var::r, side::left));
	double numx = sqrtrleft*side_val(solleft, var::vx, side::right) \
		+ sqrtrright*side_val(solright, var::vx, side::left);
	double numy = sqrtrleft*side_val(solleft, var::vy, side::right) \
		+ sqrtrright*side_val(solright, var::vy, side::left);
	double numz = sqrtrleft*side_val(solleft, var::vz, side::right) \
		+ sqrtrright*side_val(solright, var::vz, side::left);
	double den = sqrtrleft + sqrtrright;
	return{ numx/den, numy/den, numz/den };
	//return 0.125*(sqr(side_val(Uright, Vright, var::vx, side::left) + side_val(Uleft, Vleft, var::vx, side::right)) + \
				//sqr(side_val(Uright, Vright, var::vy, side::left) + side_val(Uleft, Vleft, var::vy, side::right)) + \
				//sqr(side_val(Uright, Vright, var::vz, side::left) + side_val(Uleft, Vleft, var::vz, side::right)));
}
//скорость на €чейке
vector<double> ProblemMHD1D::v(const vector<vector<double>>& sol) const
{
	double numx = val(sol, var::vx);
	double numy = val(sol, var::vy);
	double numz = val(sol, var::vz);
	return{ numx , numy , numz };
	//return 0.125*(sqr(side_val(Uright, Vright, var::vx, side::left) + side_val(Uleft, Vleft, var::vx, side::right)) + \
					//sqr(side_val(Uright, Vright, var::vy, side::left) + side_val(Uleft, Vleft, var::vy, side::right)) + \
					//sqr(side_val(Uright, Vright, var::vz, side::left) + side_val(Uleft, Vleft, var::vz, side::right)));
}



//квадрат скорости на границе (left)-й и (right)-й €чеек
double ProblemMHD1D::v2_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const
{
	/*
	double sqrtrleft = sqrt(side_val(Uleft, Vleft, var::r, side::right));
	double sqrtrright = sqrt(side_val(Uright, Vright, var::r, side::left));
	double numx = sqrtrleft*side_val(Uleft, Vleft, var::vx, side::right) \
		+ sqrtrright*side_val(Uright, Vright, var::vx, side::left);
	double numy = sqrtrleft*side_val(Uleft, Vleft, var::vy, side::right) \
		+ sqrtrright*side_val(Uright, Vright, var::vy, side::left);
	double numz = sqrtrleft*side_val(Uleft, Vleft, var::vz, side::right) \
		+ sqrtrright*side_val(Uright, Vright, var::vz, side::left);
	double den = sqrtrleft + sqrtrright;
	return (sqr(numx) + sqr(numy) + sqr(numz)) / sqr(den);
	*/

	vector<double> vav = v_av(solleft, solright);
	return vav[0] * vav[0] + vav[1] * vav[1] + vav[2] * vav[2];

	//return 0.125*(sqr(side_val(Uright, Vright, var::vx, side::left) + side_val(Uleft, Vleft, var::vx, side::right)) + \
			//sqr(side_val(Uright, Vright, var::vy, side::left) + side_val(Uleft, Vleft, var::vy, side::right)) + \
			//sqr(side_val(Uright, Vright, var::vz, side::left) + side_val(Uleft, Vleft, var::vz, side::right)));
}
//квадрат скорости на €чейке
double ProblemMHD1D::v2(const vector<vector<double>>& sol) const
{
	vector<double> vav = v(sol);
	return vav[0] * vav[0] + vav[1] * vav[1] + vav[2] * vav[2];

}

//”точненна€ скорость звука на разрыве 
double ProblemMHD1D::d_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const
{
	double rleft = side_val(solleft, var::r, side::right);
	double rright = side_val(solright, var::r, side::left);
	double sqrtrleft = sqrt(rleft);
	double sqrtrright = sqrt(rright);
	
	double uleft = side_val(solleft, var::vx, side::right);
	double uright = side_val(solright, var::vx, side::left);
	double pleft = side_val(solleft, var::p, side::right);
	double pright = side_val(solright, var::p, side::left);
    double c2left = (gamma * pleft)/rleft;
    double c2right = (gamma * pright)/rright;

	double w = 0.5*(sqrtrleft*sqrtrright) / sqr(sqrtrleft + sqrtrright);
		
    return sqrt((sqrtrleft*c2left + sqrtrright*c2right) / (sqrtrleft + sqrtrright) + w*sqr(uright-uleft));
}

//напр€жЄнности ћѕ на границе €чеек
vector<double> ProblemMHD1D::H_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const
{
	return{ Hx, 0.5*(side_val(solright, var::Hy, side::left) + side_val(solleft, var::Hy, side::right)), 0.5*(side_val(solright, var::Hz, side::left) + side_val(solleft, var::Hz, side::right)) };
	
}
//напр€жЄнности ћѕ на €чейке
vector<double> ProblemMHD1D::H(const vector<vector<double>>& sol) const
{
	return{ Hx, val(sol, var::Hy), val(sol, var::Hz) };
}
//квадрат напр€жЄнности ћѕ на границе (left)-й и (right)-й €чеек
double ProblemMHD1D::H2_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const
{

	vector<double> Hav = H_av(solleft, solright);
	return Hav[0] * Hav[0] + Hav[1] * Hav[1] + Hav[2] * Hav[2];
}
//квадрат напр€жЄнности ћѕ на €чейке
double ProblemMHD1D::H2(const vector<vector<double>>& sol) const
{
	vector<double> Hav = H(sol);
	return Hav[0] * Hav[0] + Hav[1] * Hav[1] + Hav[2] * Hav[2];

}
// ѕолное давление
double ProblemMHD1D::PT(const vector<vector<double>>& sol) const
{
	return val(sol, var::p) + H2(sol) / 2;
}


// Ќекоторые коэффициенты дл€ матриц
double ProblemMHD1D::B2(const vector<vector<double>>& sol) const
{
	double a = val(sol, var::Hy) * val(sol, var::Hy) + val(sol, var::Hz) * val(sol, var::Hz);
	if (a == 0)
		return 1 / sqrt(2);
	else
		return val(sol, var::Hy) / sqrt(a);
}
double ProblemMHD1D::B3(const vector<vector<double>>& sol) const
{
	double a = val(sol, var::Hy) * val(sol, var::Hy) + val(sol, var::Hz) * val(sol, var::Hz);
	if (a == 0)
		return 1 / sqrt(2);
	else
		return val(sol, var::Hz) / sqrt(a);
}
double ProblemMHD1D::S1(const vector<vector<double>>& sol) const
{
	if (Hx >= 0) return 1;
	else return -1;
}
double ProblemMHD1D::AlfaF(const vector<vector<double>>& sol) const
{
	double cc = c(sol);
	double cs = As(sol);
	double cf = Af(sol);

	if (cs == cf) return 1 / sqrt(2);
	else
		return (sqrt(cc*cc - cs*cs)) / (sqrt(cf*cf - cs*cs));
}
double ProblemMHD1D::AlfaS(const vector<vector<double>>& sol) const
{
	double cc = c(sol);
	double cs = As(sol);
	double cf = Af(sol);

	if (cs == cf) return 1 / sqrt(2);
	else
		return (sqrt(cf*cf - cc*cc)) / (sqrt(cf*cf - cs*cs));
}

// —корости распространени€

// Ќаибольша€ скорость распространени€ на временном уровне
double ProblemMHD1D::CFLSpeedMax(const vector<vector<vector<double>>>& SOL) const
{
	double umax = 0.0;  // variable for the max value of normal speed v_x = u
	double AfMax = 0.0; // variable for the max value of the fast magnetosonic speed a_f
	double u = 0.0, A_f = 0.0; // buffers

	const int nx = ptrprm->nx;

	// Searching for the corresponding max speeds in elements
	for (int cell = 0; cell < nx; ++cell) // до nx+1???
	{
		u = fabs(val(SOL[cell], var::vx));
		A_f = fabs(Af(SOL[cell]));

		if (u > umax) umax = u;
		if (A_f > AfMax) AfMax = A_f;
	}// for i

	// Computation of the max velocity of discontinuities
	return umax + AfMax;
}
// јльфвеновка€ скорость
double ProblemMHD1D::SpeedAlfven(const vector<vector<double>>& sol) const
{
	return fabs(Hx) / sqrt(val(sol, var::r));
}
// Ѕыстра€ магнитозвукова€ скорость															
double ProblemMHD1D::Af(const vector<vector<double>>& sol) const
{
	double cc = c(sol)*c(sol);
	double aa = SpeedAlfven(sol)*SpeedAlfven(sol);
	double HH = H2(sol);
	double r = val(sol, var::r); 

	return sqrt(0.5*cc + 0.5* HH / r + 0.5*sqrt((cc + HH / r) * (cc + HH / r) - 4 * cc * aa));
}
// ћедленна€ магнитозвукова€ скорость
double ProblemMHD1D::As(const vector<vector<double>>& sol) const
{
	double cc = c(sol)*c(sol);
	double aa = SpeedAlfven(sol)*SpeedAlfven(sol);
	double HH = H2(sol);
	double r = val(sol, var::r);

	return sqrt(0.5*cc + 0.5* HH / r - 0.5*sqrt((cc + HH / r) * (cc + HH / r) - 4 * cc * aa));
}

//«начени€ левых и правых собственных векторов на всей сетке
////(на левых границах €чеек - между (i-1)-й и i-й €чейками)
//void ProblemGas1D::omega(const vector<vector<vector<double>>>& SOL, \
//	vector<vector<vector<double>>>& LW, \
//	vector<vector<vector<double>>>& RW, \
//	const initializer_list<int>& list) const
//{
//	int nx = ptrprm->nx;
//
//	//ƒл€ временного хранени€ собств. векторов
//	vector<double> rw, lw;
//
//	for (int cell = 0; cell < nx + 1; ++cell)
//	{
//		vector<vector<double>>& rw = RW[cell];
//		vector<vector<double>>& lw = LW[cell];
//
//		const vector<vector<double>>& solright = cell < nx ? SOL[cell]     : SOL[nx + 1];
//		const vector<vector<double>>& solleft  = cell > 0  ? SOL[cell - 1] : SOL[nx];
//		
//		double cav = c_av(solleft, solright);
//		double hav = h_av(solleft, solright);
//		double v2av = v2_av(solleft, solright);
//
//		double vxr = side_val(solright, var::vx, side::left);
//		double vyr = side_val(solright, var::vy, side::left);
//		double vzr = side_val(solright, var::vz, side::left);
//
//		double vxl = side_val(solleft, var::vx, side::right);
//		double vyl = side_val(solleft, var::vy, side::right);
//		double vzl = side_val(solleft, var::vz, side::right);
//
//		const double gm = gamma - 1.0;
//
//		rw[0] = { 1, 0, 0, 1, 1 };
//		rw[1] = { -cav + 0.5*(vxl + vxr), 0, 0, 0.5*(vxl + vxr), cav + 0.5*(vxl + vxr) };
//		rw[2] = { 0.5*(vyl + vyr), 1, 0, 0.5*(vyl + vyr), 0.5*(vyl + vyr) };
//		rw[3] = { 0.5*(vzl + vzr), 0, 1, 0.5*(vzl + vzr), 0.5*(vzl + vzr) };
//		rw[4] = { hav - 0.5*cav*(vxl + vxr), 0.5*(vyl + vyr), 0.5*(vzl + vzr), -sqr(cav) / gm + hav, hav + 0.5*cav*(vxl + vxr) };
//
//
//		lw[0] = { 0.5*cav*(vxl + vxr) / gm + v2av, -cav / gm - 0.5*(vxl + vxr), -0.5*(vyl + vyr), -0.5*(vzl + vzr), 1 };
//		lw[1] = { -sqr(cav)*(vyl + vyr) / gm, 0, 2.0*sqr(cav) / gm, 0, 0 };
//		lw[2] = { -sqr(cav)*(vzl + vzr) / gm, 0, 0, 2.0*sqr(cav) / gm, 0 };
//		lw[3] = { 2.0*hav - 4.0*v2av, vxl + vxr, vyl + vyr, vzl + vzr, -2 };
//		lw[4] = { -0.5*cav*(vxl + vxr) / gm + v2av, cav / gm - 0.5*(vxl + vxr), -0.5*(vyl + vyr), -0.5*(vzl + vzr), 1 };
//
//		double cft = 0.5*gm / sqr(cav);
//
//		for (int i = 0; i < 5; ++i)
//		for (int j = 0; j < 5; ++j)
//			lw[i][j] *= cft;
//
//		for (size_t j = 0; j < list.size(); ++j)
//		{
//			LW[cell][j] = lw[*(list.begin() + j)];
//			RW[cell][j] = rw[*(list.begin() + j)];
//		}
//	}//for cell
//}
// ћатрицы собственных векторов на €чейке.
void ProblemMHD1D::EigenMatricies(const vector<vector<double>>& sol, \
	vector<vector<double>>& LL, \
	vector<vector<double>>& RR) const
{

	double r  = val(sol, var::r);
	double cc = c(sol);
	double c2 = cc*cc;
	double cf = Af(sol);
	double cs = As(sol);
	double af = AlfaF(sol);
	double as = AlfaS(sol);
	double b2 = B2(sol);
	double b3 = B3(sol);
	double s1 = S1(sol);
	double d1 = sqrt(0.5*r); // = d*sqrt(2), d=sqrt(r)/2.0
	double d3 = sqrt(2.0*r); // = 2*sqrt(2)*d, d=sqrt(r)/2.0
	double d2 = 2 * sqrt(r); // = 4*d
	double sqrtr = sqrt(r);
	double sqrt2 = sqrt(2);

	RR[0] = { 1, 0, 0, r*af, r*af, r*as, r*as };
	RR[1] = { 0, 0, 0, af*cf, -af*cf, as*cs, -as*cs };
	RR[2] = { 0, -b3 / sqrt2, -b3 / sqrt2, -as*cs*b2*s1, as*cs*b2*s1, af*cf*b2*s1, -af*cf*b2*s1 };
	RR[3] = { 0, b2 / sqrt2, b2 / sqrt2, -as*cs*b3*s1, as*cs*b3*s1, af*cf*b3*s1, -af*cf*b3*s1 };
	RR[4] = { 0, 0, 0, r*af*c2, r*af*c2, r*as*c2, r*as*c2 };
	RR[5] = { 0, b3*s1*d1, -b3*s1*d1, as*cc*b2*sqrtr, as*cc*b2*sqrtr, -af*b2*cc*sqrtr, -af*b2*cc*sqrtr };
	RR[6] = { 0, -b2*s1*d1, b2*s1*d1, as*b3*cc*sqrtr, as*b3*cc*sqrtr, -af*b3*cc*sqrtr, -af*b3*cc*sqrtr };

	LL[0] = { 1, 0, 0, 0, -1 / (c2), 0, 0 };
	LL[1] = { 0, 0, -b3 / sqrt2, b2 / sqrt2, 0, b3*s1 / d3, -b2*s1 / d3 };
	LL[2] = { 0, 0, -b3 / sqrt2, b2 / sqrt2, 0, -b3*s1 / d3, b2*s1 / d3 };
	LL[3] = { 0, af*cf / (2*c2), -as*cs*b2*s1 / (2*c2), -as*cs*b3*s1 / (2*c2), af / (2*c2*p), as*b2 / (cc*d2), as*b3 / (cc*d2) };
	LL[4] = { 0, -af*cf / (2*c2), as*cs*b2*s1 / (2*c2), as*cs*b3*s1 / (2*c2), af / (2*c2*p), as*b2 / (cc*d2), as*b3 / (cc*d2) };
	LL[5] = { 0, as*cs / (2*c2), af*cf*b2*s1 / (2*c2), af*cf*b3*s1 / (2*c2), as / (2*c2*p), -af*b2 / (cc*d2), -af*b3 / (cc*d2) };
	LL[6] = { 0, -as*cs / (2*c2), -af*cf*b2*s1 / (2*c2), -af*cf*b3*s1 / (2*c2), as / (2*c2*p), -af*b2 / (cc*d2), -af*b3 / (cc*d2) };
}


//«начени€ левых и правых собственных чисел на всей сетке 
//(на левых границах €чеек - между (i-1)-й и i-й €чейками)
void ProblemMHD1D::lambda(const vector<vector<vector<double>>>& SOL, \
	const SoundVelType soundveltype, \
	vector<vector<double>>& LL, \
	const initializer_list<int>& list) const
{
	//ƒл€ временного хранени€ собств. чисел
	vector<double> lam;

	const int nx = ptrprm->nx;

	for (int cell = 0; cell < nx + 1; ++cell)
	{
		//const vector<vector<double>>& solright = cell < nx ? SOL[cell]     : SOL[nx + 1];
		//const vector<vector<double>>& solleft  = cell > 0  ? SOL[cell - 1] : SOL[nx];
		//		
		//double cav = 0;

		//switch (soundveltype)
		//{
		//case RoeSoundVel:
		//	cav = c_av(solleft, solright);
		//	break;
		//case EinfeldtSoundVel:
		//	cav = d_av(solleft, solright);
		//	break;
		//default:
		//	cav = c_av(solleft, solright);
		//}		

		////double vav = sqrt(v2_av(Uleft, Vleft, Uright, Vright));
		////lam = { vav - cav, vav, vav, vav, vav + cav };
		//
		//double vav = v_av(solleft, solright)[0];
		//		
		////if (vav >= 0)
		//	lam = { vav - cav, vav, vav, vav, vav + cav };
		////else 
		////	lam = { vav + cav, vav, vav, vav, vav - cav };

		//for (size_t j = 0; j < list.size(); ++j)
		//	LL[cell][j] = lam[*(list.begin()+j)];	

		
		double vv = v(SOL[cell])[0];
		double cs = As(SOL[cell]);
		double ca = SpeedAlfven(SOL[cell]);
		double cf = Af(SOL[cell]);
				
		//if (vav >= 0)
		lam = { vv - cf, vv - ca, vv - cs, vv, vv + cs, vv + ca, vv + cf };
		//else 
		//	lam = { vav + cav, vav, vav, vav, vav - cav };

		for (size_t j = 0; j < list.size(); ++j)
			LL[cell][j] = lam[*(list.begin()+j)];

	}//for cell
}


//Ќачальные услови€
vector<double> ProblemMHD1D::initial_var(const double x) const
{	
	double rh = init[0](x);
	double vx = init[1](x);
	double vy = init[2](x);
	double vz = init[3](x);
	double Hy = init[5](x);
	double Hz = init[6](x);
	return{ rh, rh*vx, rh*vy, rh*vz, \
		init[4](x) / (gamma - 1) + 0.5*rh*(sqr(vx) + sqr(vy) + sqr(vz)) + 0.5*(sqr(Hx) + sqr(Hy) + sqr(Hz)), \
		Hy, Hz };
}//initial_var

void ProblemMHD1D::more_to_file(ostream& str, \
	const vector<vector<vector<double>>>& SOL, int cell) const
{
	double p_i = (gamma - 1) * (SOL[cell][0][var::e] - \
		0.5*(sqr(SOL[cell][0][var::rvx]) + \
		sqr(SOL[cell][0][var::rvy]) + \
		sqr(SOL[cell][0][var::rvz])) / SOL[cell][0][var::r] - \
		0.5*(Hx*Hx + sqr(SOL[cell][0][var::Hy]) + sqr(SOL[cell][0][var::Hz])));

	str << p_i / (gamma - 1) << " " << 0.0 << " " << 0.0 << " ";

}
