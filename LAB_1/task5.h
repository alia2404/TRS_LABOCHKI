#pragma once
#include <string>
void t5(int n, double x_start, double x_fin, std::string out_filename, double& x_step);
void shotgun_method(int n, double x_start, double x_fin, std::string out_filename);
double dudt(double x, double u, double v);
double dvdt(double x, double u, double v);
double rk4_search(int n, double x_start, double x_fin, double u_start, std::string out_filename);
double psi(double v);
