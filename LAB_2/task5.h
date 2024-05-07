#pragma once
#include <string>
#include "my_data.h"

void t5(double x_step, double t_step, double t_fin, double delta, std::string out_filename);
void t5_u(double x_step, double t_step, double t_fin, double delta, std::string out_filename);
double calc_delta(std::vector<double> v1, std::vector<double> v2);
