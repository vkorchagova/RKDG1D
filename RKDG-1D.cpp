// RKDG-1D.cpp: определ€ет точку входа дл€ консольного приложени€.
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "defs.h"
#include "Params.h"
#include "Mesh1D.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
	//»нициализируем параметры задачи
	//Params prm(LeftTriangle);
	Params prm(Sod);
	//Params prm(BrioWu);
	//Params prm(Lax);
	//Params prm(Rarefaction);
	//Params prm(leftWoodward);
	//Params prm(rightWoodward);
    //Params prm(Shock);
	//Params prm(TestSmallTravelingWave);
	
	//»нициализируем сетку и сохран€ем решение в начальный момент времени
	Mesh1D mesh(prm);

	//—обственно цикл по времени
	do
	{
		mesh.TimeStepping();
	}
	while (!mesh.finished());
	
    return 0;
}
