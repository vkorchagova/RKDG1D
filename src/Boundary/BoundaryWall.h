#ifndef BOUNDARYWALL_H_
#define BOUNDARYWALL_H_

#include <vector>

#include "Boundary.h"

#define Wall BoundaryWall 

class BoundaryWall :
	public Boundary
{
protected:

public:
	BoundaryWall(const BaseParams& prm, const Problem& prb);
	~BoundaryWall();

	//–≈јЋ»«ј÷»я ¬»–“”јЋ№Ќќ… ‘”Ќ ÷»» 
	//—обственно, происходит формирование √” путем создани€ решени€ в фиктивной €чейке
	//заполнение значений новых решений U и наклонов V
	//на левой (nx) и правой (nx+1) фиктивных €чейках
	//рассчитываетс€ по решению (UU) и наклонам (VV) на всех €чейках сетки
	void FillVirtualCells(const vector<vector<vector<double>>>& SOL) const;
};

#endif