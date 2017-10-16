#ifndef BOUNDARYSOFT_H_
#define BOUNDARYSOFT_H_

#include <vector>

#include "Boundary.h"

#define Soft BoundarySoft 

class BoundarySoft :
	public Boundary
{
protected:

public:
	BoundarySoft(const BaseParams& prm, const Problem& prb);
	~BoundarySoft();

	//–≈јЋ»«ј÷»я ¬»–“”јЋ№Ќќ… ‘”Ќ ÷»» 
	//—обственно, происходит формирование √” путем создани€ решени€ в фиктивной €чейке
	//заполнение значений новых решений U и наклонов V
	//на левой (nx) и правой (nx+1) фиктивных €чейках
	//рассчитываетс€ по решению (UU) и наклонам (VV) на всех €чейках сетки
	void FillVirtualCells(const vector<vector<vector<double>>>& SOL) const;
};

#endif