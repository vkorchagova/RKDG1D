// RKDG-1D.cpp: определяет точку входа для консольного приложения.
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
    //Инициализируем параметры задачи
    //Params prm(LeftTriangle);
    Params prm(Sod);
    //Params prm(BrioWu);
    //Params prm(Lax);
    //Params prm(Rarefaction);
    //Params prm(leftWoodward);
    //Params prm(rightWoodward);
    //Params prm(Shock);
    //Params prm(TestSmallTravelingWave);

    //Инициализируем сетку и сохраняем решение в начальный момент времени
    Mesh1D mesh(prm);

    //Собственно цикл по времени
    do
    {
        mesh.TimeStepping();
    }
    while (!mesh.finished());

    return 0;
}

