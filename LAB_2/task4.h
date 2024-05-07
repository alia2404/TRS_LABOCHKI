#pragma once
#include <string>
#include "my_data.h"

/*
explicit method
*/
void t4(double x_step, double t_step, double t_fin, double sigma, std::string out_filename);
void explicit_t4(double x_step, double t_step, double t_fin, std::string out_filename);
void implicit_t4(double x_step, double t_step, double t_fin, double sigma, std::string out_filename);
