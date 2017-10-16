#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include <vector>

#include "defs.h"
#include "Problem.h"
#include "ProblemGas1D.h"
#include "ProblemMHD1D.h"
#include "ProblemTransfer1D.h"

using namespace std;

#define SET_NO_CONST_REF(u, constu) vector<vector<double>>& u = const_cast<vector<vector<double>>&>(constu);

class Boundary
{
protected:
	//”казатель на объект, содержащий параметры задачи
	const BaseParams* ptrprm;

	//”казатель на объект, содержащий решаемую задачу
	const Problem* ptrprb;

	// опирование решени€ (U, V, W) -> (copyU, copyV, copyW)
	void CopySol(const vector<vector<double>>& SOL, vector<vector<double>>& copySOL) const;

public:
	//¬»–“”јЋ№Ќјя ‘”Ќ ÷»я 
	//—обственно, происходит формирование √” путем создани€ решени€ в фиктивной €чейке
	//заполнение значений новых решений U и наклонов V
	//на левой (nx) и правой (nx+1) фиктивных €чейках
	//рассчитываетс€ по решению (UU) и наклонам (VV) на всех €чейках сетки
	virtual void FillVirtualCells(const vector<vector<vector<double>>>& SOL) const = 0;

	// онструктор (prm - объект, содержащий параметры задачи
	//prb - решаема€ задача
	Boundary(const BaseParams& prm, const Problem& prb);

	//ƒеструктор
	~Boundary();
};

#endif
