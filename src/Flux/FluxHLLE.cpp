#include "FluxHLLE.h"

//КОНСТРУКТОР HLLE
FluxHLLE::FluxHLLE(const BaseParams& prm, Problem& prb) : Flux(prm, prb)
{
    int nx = ptrprm->nx;
    
	//Заготовки под L
    L.resize(nx + 1);
    for (int row = 0; row < nx + 1; ++row)
    {
        L[row].resize(2);
    }
    
	//Заготовки под необходимые (временные) векторы
    LdU.resize(5), RdU.resize(5), hllL.resize(5), hllR.resize(5);
}

//ДЕСТРУКТОР HLLE
FluxHLLE::~FluxHLLE()
{};

void FluxHLLE::step(const vector<vector<double>>& UU, const vector<vector<double>>& VV, \
                   vector<vector<double>>& DU, vector<vector<double>>& DV, \
                   const double cft)
{
    int nx = ptrprm->nx;
    double h = ptrprm->h;
    int dim = ptrprb->dim;
    
	//Находим конвективные потоки
    ptrprb->convFlux(UU, VV);
    
	//Находим собственные числа
	ptrprb_toGas->lambda(UU, VV, EinfeldtSoundVel, L, { 0, dim - 1 });
    
	//Основной цикл по ячейкам
    for (int cell = 0; cell < nx; ++cell)
    {
		//Установка ссылок на решения на трех ячейках	:
		//  myu(), myv() - на своей ячейке
		//  leftu(), leftv() - на соседней слева ячейке
		//  rightu(), rightv() - на соседней справа ячейке 
		//и на конв.потоки на трех ячейках:
		//  myflux(), myfluxL(), myfluxR() - конв.потоки на своей ячейке (по центру, слева, справа)
		//  rightfluxL(), leftfluxR() - конв.потоки на правой и левой соседних ячейках (слева, справа)
        setlocsolflux(UU, VV, cell);
        
		//Скачок значений решения на ячейках
        for (size_t val = 0; val < 5; ++val)
        {
            LdU[val] = side_val(myu(), myv(), (var)val, side::left)       - side_val(leftu(), leftv(), (var)val, side::right);
            RdU[val] = side_val(rightu(), rightv(), (var)val, side::left) - side_val(myu(), myv(), (var)val, side::right);
        }
        
		
		//Проверка направлений переноса
		if (L[cell][0] > 0)
			for (size_t val = 0; val < 5; ++val)
				hllL[val] = leftfluxR()[val];
		else
			if (L[cell][1] < 0)
				for (size_t val = 0; val < 5; ++val)
					hllL[val] = myfluxL()[val];
		else
			for (size_t val = 0; val < 5; ++val)
				hllL[val] = (L[cell][1] * leftfluxR()[val] - L[cell][0] * myfluxL()[val] + L[cell][1] * L[cell][0] * LdU[val]) / (L[cell][1] - L[cell][0]);

		//Проверка направлений переноса
		if (L[cell + 1][0] > 0)
			for (size_t val = 0; val < 5; ++val)
				hllR[val] = myfluxR()[val];
		else
			if (L[cell + 1][1] < 0)
				for (size_t val = 0; val < 5; ++val)
					hllR[val] = rightfluxL()[val];
		else
			for (size_t val = 0; val < 5; ++val)
				hllR[val] = (L[cell + 1][1] * myfluxR()[val] - L[cell + 1][0] * rightfluxL()[val] + L[cell + 1][1] * L[cell + 1][0] * RdU[val]) / (L[cell + 1][1] - L[cell + 1][0]);



        for (size_t val = 0; val < 5; ++val)
        {
            hllL[val] = (L[cell][1] * leftfluxR()[val] - L[cell][0] * myfluxL()[val]      + L[cell][1] * L[cell][0] * LdU[val]) / (L[cell][1] - L[cell][0]);
            hllR[val] = (L[cell+1][1] * myfluxR()[val] - L[cell+1][0] * rightfluxL()[val] + L[cell+1][1] * L[cell+1][0] * RdU[val]) / (L[cell+1][1] - L[cell+1][0]);
        }
        
		//Пересчитываем средние значения и потоки
        for (size_t val = 0; val < 5; ++val)
        {
            DU[cell][val] = -(cft / h) * (hllR[val] - hllL[val]);
            DV[cell][val] = -3.0*(cft / h) * (-2.0*(myfluxL()[val] + 4.0*myflux()[val] + myfluxR()[val]) / 6.0 + hllR[val] + hllL[val]);
        }
    }
}

