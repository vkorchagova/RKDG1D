#ifndef DEFS_H_
#define DEFS_H_

#include <vector>
#include <iostream>
#include <cstring>
#include <cmath>
#include <functional>

using namespace std;

//Тип данных "рабочая переменная"
enum var { r = 0, U = 0, \
	       rvx = 1, rvy = 2, rvz = 3, e = 4, \
		   Hy = 5, Hz = 6, \
		   vx = 101, vy = 102, vz = 103, p = 104 };

//Способы вычисления скорости звука
enum SoundVelType { RoeSoundVel = 0, EinfeldtSoundVel = 1, SemisumSoundVel = 2 };

typedef struct
{
	double L;     // Длина области течения
	double T;     // Расчетное время моделирования

	int nx;       // Число шагов по пространству
	double Co;    // Число Куранта
	int deltacnt; // Шаг вывода в файл

	double h;     // Шаг по пространству
	double tau;   // Шаг по времени
} BaseParams;

//Тип данных "граничное значение"
enum side { left, right };

//Возведение числа в квадрат
inline double sqr(const double p) { return p*p; };

//Прибавление к одному тройному массиву другого тройного массива
vector<vector<vector<double>>>& operator += (vector<vector<vector<double>>>& a, const vector<vector<vector<double>>>& b);

//Прибавление к одному двумерному массиву другого двумерного массива
vector<vector<double>>& operator += (vector<vector<double>>& a, const vector<vector<double>>& b);

//Вычитание из одного двумерного массива другого двумерного массива
vector<vector<double>>& operator -= (vector<vector<double>>& a, const vector<vector<double>>& b);

//Домножение трехмерного массива на число
vector<vector<vector<double>>>& operator *= (vector<vector<vector<double>>>& a, const double b);
vector<vector<vector<double>>> operator * (const vector<vector<vector<double>>>& a, const double b);

//Домножение двумерного массива на число
vector<vector<double>>& operator *= (vector<vector<double>>& a, const double b);
vector<vector<double>> operator * (const vector<vector<double>>& a, const double b);

//Домножение вектора на число
vector<double>& operator *= (vector<double>& a, const double b);
vector<double> operator * (const vector<double>& a, const double b);
vector<double> operator / (const vector<double>& a, const double b);

//Прибавление к одному вектору второго
vector<double>& operator += (vector<double>& a, const vector<double>& b);

//Вычитание из одного вектора второго
vector<double>& operator -= (vector<double>& a, const vector<double>& b);


//Умножение матрицы на вектор
void prodMatrVec(const vector<vector<double>>& A, \
    const vector<double>& b, \
    vector<double>& c);

//Умножение матриц из собственных векторов и собств. чисел (для КИР)
void prodWrAbsLWl(const vector<vector<double>>& Wr, \
    const vector<vector<double>>& Wl, \
    const vector<double>& L, \
    vector<vector<double>>& Prod);

#endif
