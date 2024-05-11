#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <stdlib.h> 

#define ITER_LIM 100

using namespace std;

vector<vector<double>> random_matrix(int M, int N);

vector<vector<double>> matrix_prod(vector<vector<double>> A, double C);

vector<vector<double>> matrix_prod(vector<vector<double>> A, vector<vector<double>> B);

vector<vector<double>> matrix_add(vector<vector<double>> A, vector<vector<double>> B);

vector<vector<double>> matrix_sub(vector<vector<double>> A, vector<vector<double>> B);

vector<vector<double>> low_tr_matrix_inv(vector<vector<double>> A);

vector<vector<double>> up_tr_matrix_inv(vector<vector<double>> A);

void back_substitution(vector<vector<double>> A, vector<double>& X, vector<double> B);

void separate_matrix(vector<vector<double>> A, vector<vector<double>>& D, vector<vector<double>>& E, vector<vector<double>>& F);

double scal_prod(vector<double> x, vector<double> y);

double eucl_norm(vector<double> x);

void matrix_transpose(vector<vector<double>> A, vector<vector<double>>& A_T);

void forward_substitution(vector<vector<double>> A, vector<double>& X, vector<double> B);

vector<vector<double>> sor_step(vector<vector<double>> X, vector<vector<double>> D, vector<vector<double>> F, vector<vector<double>> B, vector<vector<double>> DmwE_inv, double w);

vector<vector<double>> solve_sor(vector<vector<double>> A, vector<vector<double>> B, double err, double w, int& iter_number);

vector<vector<double>> solve_sor(vector<vector<double>> A, vector<vector<double>> B, double err, double w, vector<vector<double>> X_true, int& iter_number);

vector<vector<double>> sym_sor_step(vector<vector<double>> X, vector<vector<double>> D, vector<vector<double>> E, vector<vector<double>> F, vector<vector<double>> B, vector<vector<double>> DmwE_inv, vector<vector<double>> DmwF_inv, double w);

vector<vector<double>> solve_sym_sor(vector<vector<double>> A, vector<vector<double>> B, double& err, double w, int& iter_number, int mode);

vector<vector<double>> solve_sym_sor(vector<vector<double>> A, vector<vector<double>> B, double err, double w, vector<vector<double>> X_true, int& iter_number);

void s_matrix_define(vector<vector<double>> A, vector<vector<double>>& S, vector<vector<double>>& SIGN);

void Square_root_method(vector<vector<double>> A, vector<double> B, vector<double>& X);

void lu_matrix(vector<vector<double>> A, vector<vector<double>>& L, vector<vector<double>>& U);

void LU_decay(vector<vector<double>> A, vector<double> B, vector<double>& X);
