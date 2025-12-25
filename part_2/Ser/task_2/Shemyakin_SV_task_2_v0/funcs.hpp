#ifndef FUNCS_HPP
#define FUNCS_HPP

#include<cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

// Функция для вычисления теоретического решения
double ans(double A, double y0, double x);

// Схема 1: Явная схема Эйлера
double scheme_1(double y0, double A, int N);

// Схема 2: Неявная схема Эйлера
double scheme_2(double y0, double A, int N);

// Схема 3: Схема Кранка-Николсона (симметричная)
double scheme_3(double y0, double A, int N);

double scheme_4(double y0, double A, int N);

#endif 