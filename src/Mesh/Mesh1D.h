#ifndef MESH1D_H_
#define MESH1D_H_

#include <vector>
#include <sstream>
#include <fstream>

#include "defs.h"
#include "Params.h"

#include "Flux.h"
#include "Limiter.h"
#include "Timestep.h"



using namespace std;

class Mesh1D
{
private:
    const Params* ptrprm;

	//Номер текущего шага по времени
	int step_num;

	//U1, V1 --- решение и наклон на текущем/предыдущем шаге
	//U2, V2 --- решение и наклон на текущем/предыдущем шаге
	//В позициях (n) и (n+1) - левая и правая фиктивные ячейки соответственно
	vector<vector<vector<double>>> SOL1;
	vector<vector<vector<double>>> SOL2;
	  
    //ptrU,    ptrV     ---  указатели на решение и наклоны на предыдущем шаге
	//ptrUnew, ptrVnew  ---  указатели на решение и наклоны на текушем шаге
	vector<vector<vector<double>>> *ptrSOL;
	vector<vector<vector<double>>> *ptrSOLnew;
	
	//Поток вывода, связанный с текстовым файлом (телефайлом)
	ofstream telefile;
    
	//Инициализация телефайла
    void init_telefile();


	//Заполнение вектора решения на 0-м шаге начальными условиями
	void set_initial();

public:
	//Возвращает текущий номер шага по времени
	int step() const { return step_num; }
	
	//Возвращает указатель на решение на предыдущем шаге
    const vector<vector<vector<double>>>& SOL() const { return *ptrSOL; };

	//Возвращает указатель на решение на текущем шаге
    vector<vector<vector<double>>>& SOLnew() const { return *ptrSOLnew; };

	//Выполнение одного шага по времени
	void TimeStepping();

	//Признак того, что счет закончен
	bool finished() const { return !(step_num * ptrprm->basic.tau < ptrprm->basic.T); }
	    
	//Конструктор (prm - объект, содержащий параметры задачи)
    Mesh1D(const Params& prm_);

	//Деструктор
    ~Mesh1D();

	//Вычисление положения центра i-й ячейки сетки
    double cnt(const int i){ return (ptrprm->basic.h) * (i + 0.5); };

	//Сохранение информации о решении на текущем шаге в телефайл
    void add_to_telefile();      
};


#endif

