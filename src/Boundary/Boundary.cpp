#include "Boundary.h"

//Конструктор ГУ
Boundary::Boundary(const BaseParams& prm, const Problem& prb)
{
	ptrprm = &prm;
	ptrprb = &prb;
}

//Деструктор ГУ
Boundary::~Boundary()
{
}

//Копирование решения (U, V, W) -> (copyU, copyV, copyW)
void Boundary::CopySol(const vector<vector<double>>& SOL, vector<vector<double>>& copySOL) const
{
	copySOL = SOL;	
};

