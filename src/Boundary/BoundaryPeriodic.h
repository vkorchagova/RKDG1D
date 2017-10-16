#ifndef BOUNDARYPERIODIC_H_
#define BOUNDARYPERIODIC_H_

#include <vector>

#include "Boundary.h"

#define Periodic BoundaryPeriodic 

class BoundaryPeriodic :
	public Boundary
{
protected:

public:
	BoundaryPeriodic(const BaseParams& prm, const Problem& prb);
	~BoundaryPeriodic();

	//–≈јЋ»«ј÷»я ¬»–“”јЋ№Ќќ… ‘”Ќ ÷»» 
	//—обственно, происходит формирование √” путем создани€ решени€ в фиктивной €чейке
	//заполнение значений новых решений U и наклонов V
	//на левой (nx) и правой (nx+1) фиктивных €чейках
	//рассчитываетс€ по решению (UU) и наклонам (VV) на всех €чейках сетки
	void FillVirtualCells(const vector<vector<vector<double>>>& SOL) const;
};

#endif