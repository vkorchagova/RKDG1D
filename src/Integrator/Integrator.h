#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include <vector>
#include <functional> 

#include "defs.h"

using namespace std;

enum IntegrPoints { Gauss = 0, Lobatto = 1 };

const vector<vector<vector<double>>> pos
(
{
	//Положения гауссовых точек
	{
		{ 0.0 },
		{ -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) },
		{ -sqrt(0.6), 0.0, sqrt(0.6) },
		{ -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526 },
		{ -0.906179845938664, -0.538469310105683, 0.0, 0.538469310105683, 0.906179845938664 }
	},

	//Положения точек Гаусса-Лобатто
	{
		{ 0.0 },
		{ -1.0, 1.0 },
		{ -1.0, 0.0, 1.0 },
		{ -1.0, -1.0 / sqrt(5.0), 1.0 / sqrt(5.0), 1.0 },
		{ -1.0, -sqrt(21.0) / 7.0, 0.0, sqrt(21.0)/7.0 , 1.0}
	},
}
);

const vector<vector<vector<double>>> wgt
({
	//Веса гауссовых точек
	{
		{ 2.0 },
		{ 1.0, 1.0 },
		{ 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 },
		{ 0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538 },
		{ 0.23692688505618908, 0.47862867049936647, 0.5688888888888889, 0.47862867049936647, 0.23692688505618908 }
	},

	//Веса точек Гаусса-Лобатто
	{
		{ 2.0 },
		{ 1.0, 1.0 },
		{ 1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0 },
		{ 1.0 / 6.0, 5.0 / 6.0, 5.0 / 6.0, 1.0 / 6.0 },
		{ 0.1, 49.0/90.0, 32.0/45.0, 49.0/90.0, 0.1 }
	}
}
);

class Integrator
{
private:	
	int ngp;
	const vector<double>& weight;
	const vector<double>& position;
public:
	Integrator(IntegrPoints pts, int ngp_) : ngp(ngp_), weight(wgt[pts][ngp_ - 1]), position(pos[pts][ngp_ - 1]) {};
	virtual ~Integrator() {};

	double integrate (function<double(double)> fo)
	{
		double res = 0.0;
		for (int i = 0; i < ngp; ++i)
		{
			res += weight[i] * fo(position[i]);
		}//for i
		return 0.5*res;
	} //integrate

	vector<double> integrate(function<vector<double>(double)> fo, int dim)
	{
		vector<double> res(dim, 0.0);
		for (int i = 0; i < ngp; ++i)
		{
			res += fo(position[i]) * weight[i];
		}//for i
		res *= 0.5;
		return res;
	} //integrate
};

#endif
