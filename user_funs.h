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
matrix hessian(matrix x, matrix ud1, matrix ud2);
matrix ff4R(matrix theta, matrix ud1, matrix ud2);
matrix gradient4R(matrix theta, matrix ud1, matrix ud2);
void load_data(matrix& X, matrix& Y);
double calculate_accuracy(matrix theta, matrix X, matrix Y);



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

matrix ff5R_masa(matrix x, matrix ud1, matrix ud2);
matrix ff5R_ugiecie(matrix x, matrix ud1, matrix ud2);
matrix ff5R_naprezenie(matrix x, matrix ud1, matrix ud2) ;
matrix ff5R_base(matrix x, matrix ud1, matrix ud2);
matrix ff5R(matrix x, matrix ud1, matrix ud2);
matrix f_line_powell(matrix h, matrix xk, matrix dk);
matrix ff6T(matrix x, matrix ud1, matrix ud2);
matrix ff6R(matrix x, matrix ud1, matrix ud2);
matrix df6(double t, matrix Y, matrix ud1, matrix ud2);