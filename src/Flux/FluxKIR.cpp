#include "FluxKIR.h"

//КОНСТРУКТОР КИР
FluxKIR::FluxKIR(const BaseParams& prm, Problem& prb) : Flux(prm, prb)
{
    int nx = ptrprm->nx;

    //заготовки под RW, LW
    LW.resize(nx + 1); 	RW.resize(nx + 1);
    for (int row = 0; row < nx + 1; ++row)
    {
        LW[row].resize(5);
        RW[row].resize(5);
        for (int col = 0; col < 5; ++col)
        {
            LW[row][col].resize(5);
            RW[row][col].resize(5);
        }
    }

    //Заготовки под L
    L.resize(nx + 1);
    for (int row = 0; row < nx + 1; ++row)
    {
        L[row].resize(5);
    }

    //заготовки под прозведения RW |L| LW
    LMatr.resize(5), RMatr.resize(5);
    for (int row = 0; row < 5; ++row)
    {
        LMatr[row].resize(5);
        RMatr[row].resize(5);
    }

    //Заготовки под необходимые (временные) векторы
    LdU.resize(5), RdU.resize(5), ProdL.resize(5), ProdR.resize(5);
}//KIR::KIR

//ДЕСТРУКТОР КИР
FluxKIR::~FluxKIR()
{};


void FluxKIR::step(const vector<vector<double>>& UU, const vector<vector<double>>& VV, \
	vector<vector<double>>& DU, vector<vector<double>>& DV, \
    const double cft)
{
    int nx = ptrprm->nx;
	double h = ptrprm->h;
	int dim = ptrprb->dim;

	//Находим конвективные потоки
	ptrprb->convFlux(UU, VV);

    //Находим собственные векторы и собственные числа
	ptrprb_toGas->omega(UU, VV, LW, RW, { 0, (dim - 1) / 2 - 1, (dim - 1) / 2, (dim - 1) / 2 + 1, dim - 1 });
	ptrprb_toGas->lambda(UU, VV, RoeSoundVel, L, { 0, (dim - 1) / 2 - 1, (dim - 1) / 2, (dim - 1) / 2 + 1, dim - 1 });

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
            LdU[val] = side_val(myu(), myv(), (var)val, side::left) - side_val(leftu(), leftv(), (var)val, side::right);
            RdU[val] = side_val(rightu(), rightv(), (var)val, side::left) - side_val(myu(), myv(), (var)val, side::right);
        }

        //Перемножаем их, результат - в матрицах LMatr, RMatr
        prodWrAbsLWl(RW[cell], LW[cell], L[cell], LMatr);
        prodWrAbsLWl(RW[cell + 1], LW[cell + 1], L[cell + 1], RMatr);

        //Перемножаем на скачок решения
        prodMatrVec(LMatr, LdU, ProdL);
        prodMatrVec(RMatr, RdU, ProdR);

        //Пересчитываем средние значения и потоки
        for (size_t val = 0; val < 5; ++val)
        {
            DU[cell][val] = -0.5*(cft / h)*(((rightfluxL()[val] + myfluxR()[val]) - (myfluxL()[val] + leftfluxR()[val])) \
                - (ProdR[val] - ProdL[val]));
            
			//По формуле центральных прямоугольников
			//DV[cell][val] = -1.5*(cft / h)*(-4.0*myflux[val] + ((rightfluxL[val] + myfluxR[val]) + (myfluxL[val] + leftfluxR[val])) \
            //    - (ProdR[val] + ProdL[val]));

			//По формуле Симпсона
			DV[cell][val] = -1.5*(cft / h)*(-4.0*(myfluxL()[val]+4.0*myflux()[val]+myfluxR()[val])/6.0 + ((rightfluxL()[val] + myfluxR()[val]) + (myfluxL()[val] + leftfluxR()[val])) \
			                - (ProdR[val] + ProdL[val]));
        }

    }//for cell

}//KIR::step
