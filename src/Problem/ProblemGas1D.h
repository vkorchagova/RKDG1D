#ifndef PROBLEMGAS1D_H_
#define PROBLEMGAS1D_H_

#define Gas1D ProblemGas1D

#include "Problem.h"
class ProblemGas1D :
	public Problem
{
private:
	double gamma;

	//скорость звука на границе (left)-й и (right)-й €чеек
	double c_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const;
	double c_semisum(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const;
	//скорость звука в i-й €чейке
	double c(const vector<vector<double>>& sol) const;

	//Ёнтальпи€ звука на границе (left)-й и (right)-й €чеек
	double h_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const;
	//Ёнтальпи€ звука на i-й €чейке
	double h(const vector<vector<double>>& sol) const;

	//скорость на границе €чеек
	vector<double> v_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const;
	//скорость на i-й €чейке
	vector<double> v(const vector<vector<double>>& sol) const;

	// вадрат скорости на границе (left)-й и (right)-й €чеек
	double v2_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const;
	// вадрат скорости на i-й €чейке
	double v2(const vector<vector<double>>& sol) const;

	//”точненна€ скорость звука на разрыве
	double d_av(const vector<vector<double>>& solleft, const vector<vector<double>>& solright) const;
   
protected:

	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я
	//¬ычисление конвективного потока по заданному вектору решени€
	virtual void getFlux(const vector<double>& U, vector<double>& Flux) const;
	virtual vector<double> getFlux(const vector<double>& U) const;

public:
	ProblemGas1D(const BaseParams& prm, int order, double gam, vector<std::function<double(const double)>> initv);
	~ProblemGas1D();

	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я
	//ѕредельное значение (sd = left/right) q-й компоненты решени€,
	//рассчитываемое по решению (U) и наклону (V)
	virtual inline double side_val(const vector<vector<double>>& sol, var q, side sd) const;
	virtual inline double val(const vector<vector<double>>& sol, var q) const;
	virtual inline vector<double> gauss_val(const vector<vector<double>>& sol, double Lcoord) const;
	
	//«начени€ левых (LW) и правых (RW) собственных векторов на всей сетке 
	//(на левых границах €чеек - между (i-1)-й и i-й €чейками
	//имеющие номера, указанные в списке инициализации (по возрастанию)
	//рассчитываютс€ по решени€м (UU) и наклонам (VV) на всей сетке
	//(дл€ метода  »–)
	void omega(const vector<vector<vector<double>>>& SOL, \
		vector<vector<vector<double>>>& LW, \
		vector<vector<vector<double>>>& RW, \
		const initializer_list<int>& list) const;
	// ћатрицы собственных векторов на €чейке
	void EigenMatricies(const vector<vector<double>>& sol, \
		vector<vector<double>>& LL, \
		vector<vector<double>>& RR) const;
	
	//«начени€ собственных чисел (LL) на всей сетке  
	//имеющие номера, указанные в списке инициализации (по возрастанию)
	//(на левых границах €чеек - между (i-1)-й и i-й €чейками
	//рассчитываютс€ по решени€м (UU) и наклонам (VV) на всей сетке	
	void lambda(const vector<vector<vector<double>>>& SOL, \
		const SoundVelType soundveltype, \
		vector<vector<double>>& LL, \
		const initializer_list<int>& list) const;

	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я
	//¬озвращает вектор конс. переменных в точке x, вычисл€емый по заданным в объекте Param начальным услови€м
	virtual vector<double> initial_var(const double x) const;

	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я
	//ѕечать в телефайл
	virtual void more_to_file(ostream& str, \
		const vector<vector<vector<double>>>& SOL, int cell) const;
};

#endif
