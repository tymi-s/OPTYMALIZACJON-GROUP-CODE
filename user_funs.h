#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

matrix ff1T(matrix, matrix = NAN, matrix = NAN);

matrix ff2T(matrix, matrix = NAN, matrix = NAN);

matrix ff3T(matrix, matrix , matrix = NAN);
matrix df3R(double t, matrix Y, matrix ud1, matrix ud2);
matrix ff3R(matrix x, matrix ud1, matrix ud2);

matrix ff3R_base(matrix x, matrix ud1, matrix ud2);

matrix ff4T(matrix x, matrix ud1, matrix ud2);
matrix gradient(matrix x,matrix ud1,matrix ud2);
matrix f_line(matrix h, matrix xk, matrix dk);


matrix ff5T1_1(matrix x, matrix ud1, matrix ud2);
matrix ff5T2_1(matrix x, matrix ud1, matrix ud2);
matrix ff5T3_1(matrix x, matrix ud1, matrix ud2);

matrix ff5T1_10(matrix x, matrix ud1, matrix ud2);
matrix ff5T2_10(matrix x, matrix ud1, matrix ud2);
matrix ff5T3_10(matrix x, matrix ud1, matrix ud2);

matrix ff5T1_100(matrix x, matrix ud1, matrix ud2);
matrix ff5T2_100(matrix x, matrix ud1, matrix ud2);
matrix ff5T3_100(matrix x, matrix ud1, matrix ud2);

void setW(int i);
void makeW();
double getW();
