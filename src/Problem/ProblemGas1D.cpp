#include "ProblemGas1D.h"

const int dimGas1D = 5;

ProblemGas1D::ProblemGas1D(const BaseParams& prm, int order, double gam,\
	vector<std::function<double(const double)>> initv) : Problem(prm, dimGas1D, order, initv)
{
	gamma = gam;
}


ProblemGas1D::~ProblemGas1D()
{
}


void ProblemGas1D::getFlux(const vector<double>& U, vector<double>& Flux) const
{
	double r = U[var::r];
	double rvx = U[var::rvx];
	double rvy = U[var::rvy];
	double rvz = U[var::rvz];
	double e = U[var::e];

	double kin = 0.5*(rvx*rvx + rvy*rvy + rvz*rvz) / r;
	double press = (gamma - 1)*(e - kin);
	double vx = rvx / r;

	Flux[var::r] = rvx;
	Flux[var::rvx] = rvx*vx + press;
	Flux[var::rvy] = vx * rvy;
	Flux[var::rvz] = vx * rvz;
	Flux[var::e] = (e + press)*vx;
	
}

vector<double> ProblemGas1D::getFlux(const vector<double>& U) const
{
	vector<double> F(dim);
	
	double r = U[var::r];
	double rvx = U[var::rvx];
	double rvy = U[var::rvy];
	double rvz = U[var::rvz];
	double e = U[var::e];

	double kin = 0.5*(rvx*rvx + rvy*rvy + rvz*rvz) / r;
	double press = (gamma - 1)*(e - kin);
	double vx = rvx / r;

	F[var::r] = rvx;
	F[var::rvx] = rvx*vx + press;
	F[var::rvy] = vx * rvy;
	F[var::rvz] = vx * rvz;
	F[var::e] = (e + press)*vx;

	return F;
}


//«начени€ всех переменных на левой и правой границах €чеек
inline double ProblemGas1D::side_val(const vector<vector<double>>& sol, var q, side sd) const
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
			0.5*(sqr(side_val(sol, var::rvx, sd)) + \
			sqr(side_val(sol, var::rvy, sd)) + \
			sqr(side_val(sol, var::rvz, sd))) / side_val(sol, var::r, sd));

		break;
	}//switch
	return res;
}

//«начени€ всех переменных на €чейке
inline double ProblemGas1D::val(const vector<vector<double>>& sol, var q) const
{
	double res = 0.0;

	switch (q)
	{
	case r:
	case rvx:
	case rvy:
	case rvz:
	case e:
		res = sol[0][q];
		break;

	case vx:
	case vy:
	case vz:
		res = (sol[0][q - 100]) / val(sol, var::r);
		break;

	case p:
		res = (gamma - 1) * (val(sol, var::e) - \
			0.5*(sqr(val(sol, var::rvx)) + \
			sqr(val(sol, var::rvy)) + \
			sqr(val(sol, var::rvz))) / val(sol, var::r));
		break;
	}//switch
	return res;
}

inline vector<double> ProblemGas1D::gauss_val(const vector<vector<double>>& sol, double Lcoord) const
{	
	vector<double> res(sol[0]); 

	if (nshape > 1)
		res += sol[1] * Lcoord;

	if (nshape > 2)
		res += sol[2] * (1.5*Lcoord*Lcoord - 0.5);

	return res;	
}


//скорость звука на границе (left)-й и (right)-й €чеек
double ProblemGas1D::c_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const
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
double ProblemGas1D::c(const vector<vector<double>>& sol) const
{
	/*double hh = h(sol);
	double vv2 = v2(sol);														   „“ќ «ƒ≈—№ ѕ–ќ»—’ќƒ»“?!?!
	return sqrt((gamma - 1.0)*(hh - 0.5*vv2)); */     
	return sqrt( (val(sol, var::p) * gamma) / val(sol, var::r) );
}

double ProblemGas1D::c_semisum(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const
{
	return 0.5*(c(solleft) + c(solright));
}

//энтальни€ на границе (left)-й и (right)-й €чеек
double ProblemGas1D::h_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const
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
double ProblemGas1D::h(const vector<vector<double>>& sol) const
{
	return (val(sol, var::p) + val(sol, var::e)) / val(sol, var::r);
	
}

//скорость на границе €чеек
vector<double> ProblemGas1D::v_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const
{
	/*
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
	*/

    //return vector<double>({ side_val(solleft, var::vx, side::right) + side_val(solright, var::vx, side::left), \
		side_val(solleft, var::vy, side::right) + side_val(solright, var::vy, side::left),  \
		side_val(solleft, var::vz, side::right) + side_val(solright, var::vz, side::left) })*0.5;
    return { 0.5*(side_val(solleft, var::vx, side::right) + side_val(solright, var::vx, side::left)), \
             0.5*(side_val(solleft, var::vy, side::right) + side_val(solright, var::vy, side::left)), \
             0.5*(side_val(solleft, var::vz, side::right) + side_val(solright, var::vz, side::left)) };
	//return 0.125*(sqr(side_val(Uright, Vright, var::vx, side::left) + side_val(Uleft, Vleft, var::vx, side::right)) + \
				//sqr(side_val(Uright, Vright, var::vy, side::left) + side_val(Uleft, Vleft, var::vy, side::right)) + \
				//sqr(side_val(Uright, Vright, var::vz, side::left) + side_val(Uleft, Vleft, var::vz, side::right)));
}
//скорость на €чейке
vector<double> ProblemGas1D::v(const vector<vector<double>>& sol) const
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
double ProblemGas1D::v2_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const
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
double ProblemGas1D::v2(const vector<vector<double>>& sol) const
{
	vector<double> vav = v(sol);
	return vav[0] * vav[0] + vav[1] * vav[1] + vav[2] * vav[2];

}

//”точненна€ скорость звука на разрыве 
double ProblemGas1D::d_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const
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



//«начени€ левых и правых собственных векторов на всей сетке
//(на левых границах €чеек - между (i-1)-й и i-й €чейками)
void ProblemGas1D::omega(const vector<vector<vector<double>>>& SOL, \
	vector<vector<vector<double>>>& LW, \
	vector<vector<vector<double>>>& RW, \
	const initializer_list<int>& list) const
{
	int nx = ptrprm->nx;

	//ƒл€ временного хранени€ собств. векторов
	vector<double> rw, lw;

	for (int cell = 0; cell < nx + 1; ++cell)
	{
		vector<vector<double>>& rw = RW[cell];
		vector<vector<double>>& lw = LW[cell];

		const vector<vector<double>>& solright = cell < nx ? SOL[cell]     : SOL[nx + 1];
		const vector<vector<double>>& solleft  = cell > 0  ? SOL[cell - 1] : SOL[nx];
		
		double cav = c_av(solleft, solright);
		double hav = h_av(solleft, solright);
		double v2av = v2_av(solleft, solright);

		double vxr = side_val(solright, var::vx, side::left);
		double vyr = side_val(solright, var::vy, side::left);
		double vzr = side_val(solright, var::vz, side::left);

		double vxl = side_val(solleft, var::vx, side::right);
		double vyl = side_val(solleft, var::vy, side::right);
		double vzl = side_val(solleft, var::vz, side::right);

		const double gm = gamma - 1.0;

		rw[0] = { 1, 0, 0, 1, 1 };
		rw[1] = { -cav + 0.5*(vxl + vxr), 0, 0, 0.5*(vxl + vxr), cav + 0.5*(vxl + vxr) };
		rw[2] = { 0.5*(vyl + vyr), 1, 0, 0.5*(vyl + vyr), 0.5*(vyl + vyr) };
		rw[3] = { 0.5*(vzl + vzr), 0, 1, 0.5*(vzl + vzr), 0.5*(vzl + vzr) };
		rw[4] = { hav - 0.5*cav*(vxl + vxr), 0.5*(vyl + vyr), 0.5*(vzl + vzr), -sqr(cav) / gm + hav, hav + 0.5*cav*(vxl + vxr) };


		lw[0] = { 0.5*cav*(vxl + vxr) / gm + v2av, -cav / gm - 0.5*(vxl + vxr), -0.5*(vyl + vyr), -0.5*(vzl + vzr), 1 };
		lw[1] = { -sqr(cav)*(vyl + vyr) / gm, 0, 2.0*sqr(cav) / gm, 0, 0 };
		lw[2] = { -sqr(cav)*(vzl + vzr) / gm, 0, 0, 2.0*sqr(cav) / gm, 0 };
		lw[3] = { 2.0*hav - 4.0*v2av, vxl + vxr, vyl + vyr, vzl + vzr, -2 };
		lw[4] = { -0.5*cav*(vxl + vxr) / gm + v2av, cav / gm - 0.5*(vxl + vxr), -0.5*(vyl + vyr), -0.5*(vzl + vzr), 1 };

		double cft = 0.5*gm / sqr(cav);

		for (int i = 0; i < 5; ++i)
		for (int j = 0; j < 5; ++j)
			lw[i][j] *= cft;

		for (size_t j = 0; j < list.size(); ++j)
		{
			LW[cell][j] = lw[*(list.begin() + j)];
			RW[cell][j] = rw[*(list.begin() + j)];
		}
	}//for cell
}
// ћатрицы собственных векторов на €чейке.
void ProblemGas1D::EigenMatricies(const vector<vector<double>>& sol, \
	vector<vector<double>>& LL, \
	vector<vector<double>>& RR) const
{
	
	double cc = c(sol);
	double hh = h(sol);
	double vv2 = 0.5*v2(sol);//(sqr(Val(U, var::vx)) + sqr(Val(U, var::vy)) + sqr(Val(U, var::vz)));                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	double vx = val(sol, var::vx);
	double vy = val(sol, var::vy);
	double vz = val(sol, var::vz);

	const double gm = gamma - 1.0;

	RR[0] = { 1, 0, 0, 1, 1 };
	RR[1] = { -cc + vx, 0, 0, vx, cc + vx };
	RR[2] = { vy, 1, 0, vy, vy };
	RR[3] = { vz, 0, 1, vz, vz };
	RR[4] = { hh - cc*vx, vy, vz, -sqr(cc) / gm + hh, hh + cc*vx };

	LL[0] = { cc*vx / gm + vv2, -cc / gm - vx, -vy, -vz, 1 };
	LL[1] = { -2.0*sqr(cc)*vy / gm, 0, 2.0*sqr(cc) / gm, 0, 0 };
	LL[2] = { -2.0*sqr(cc)*vz / gm, 0, 0, 2.0*sqr(cc) / gm, 0 };
	LL[3] = { 2.0*hh - 4.0*vv2, 2.0*vx, 2.0*vy, 2.0*vz, -2 };
	LL[4] = { -cc*vx / gm + vv2, cc / gm - vx, -vy, -vz, 1 };

	double cft = 0.5*gm / sqr(cc);

	for (int i = 0; i < 5; ++i)
		for (int j = 0; j < 5; ++j)
			LL[i][j] *= cft;

}


//«начени€ левых и правых собственных чисел на всей сетке 
//(на левых границах €чеек - между (i-1)-й и i-й €чейками)
void ProblemGas1D::lambda(const vector<vector<vector<double>>>& SOL, \
	const SoundVelType soundveltype, \
	vector<vector<double>>& LL, \
	const initializer_list<int>& list) const
{
	//ƒл€ временного хранени€ собств. чисел
	vector<double> lam;

	const int nx = ptrprm->nx;

	for (int cell = 0; cell < nx + 1; ++cell)
	{
		const vector<vector<double>>& solright = cell < nx ? SOL[cell]     : SOL[nx + 1];
		const vector<vector<double>>& solleft  = cell > 0  ? SOL[cell - 1] : SOL[nx];
				
		double cav = 0;

		switch (soundveltype)
		{
		case SemisumSoundVel:
			cav = c_semisum(solleft, solright);
			break;
		case RoeSoundVel:
			cav = c_av(solleft, solright);
			break;
		case EinfeldtSoundVel:
			cav = d_av(solleft, solright);
			break;
		default:
			cav = c_av(solleft, solright);
		}		

		//double vav = sqrt(v2_av(Uleft, Vleft, Uright, Vright));
		//lam = { vav - cav, vav, vav, vav, vav + cav };
		
		double vav = v_av(solleft, solright)[0];
				
		//if (vav >= 0)
			lam = { vav - cav, vav, vav, vav, vav + cav };
		//else 
		//	lam = { vav + cav, vav, vav, vav, vav - cav };

		for (size_t j = 0; j < list.size(); ++j)
			LL[cell][j] = lam[*(list.begin()+j)];

	}//for cell
}


//Ќачальные услови€
vector<double> ProblemGas1D::initial_var(const double x) const
{	
	double rh = init[0](x);
	double vx = init[1](x);
	double vy = init[2](x);
	double vz = init[3](x);
	return{ rh, rh*vx, rh*vy, rh*vz, init[4](x) / (gamma - 1) + 0.5*rh*(sqr(vx) + sqr(vy) + sqr(vz)) };
}//initial_var

void ProblemGas1D::more_to_file(ostream& str, \
	const vector<vector<vector<double>>>& SOL, int cell) const
{
	double p_i = (gamma - 1) * (SOL[cell][0][var::e] - \
		0.5*(sqr(SOL[cell][0][var::rvx]) + \
		sqr(SOL[cell][0][var::rvy]) + \
		sqr(SOL[cell][0][var::rvz])) / SOL[cell][0][var::r]);

	str << p_i / (gamma - 1) << " " << 0.0 << " " << 0.0 << " ";

}
