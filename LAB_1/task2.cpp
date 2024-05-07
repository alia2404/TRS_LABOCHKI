#include <iostream>
#include <fstream>
#include "task2.h"

using namespace std;

void t2(double t_step, double t_fin, string out_filename) {
	int n = t_fin / t_step + 1;
	double t_step_sc = t_sc(t_step);
	double t_fin_sc = t_sc(t_fin);
	double* x_sc_arr; x_sc_arr = (double*) malloc(n * sizeof(double));
	double* v_sc_arr; v_sc_arr = (double*) malloc(n * sizeof(double));
	x_sc_arr[0] = x_sc(x0);
	v_sc_arr[0] = v_sc(v0);
	x_sc_arr[1] = x_sc_arr[0] + t_step_sc * dxdt_sc(x_sc_arr[0], v_sc_arr[0]);
	v_sc_arr[1] = v_sc_arr[0] + t_step_sc * dvdt_sc(x_sc_arr[0], v_sc_arr[0]);
	for (int i = 2; i < n; i++) {
		x_sc_arr[i] = x_sc_arr[i - 1]
			+ t_step_sc * (1.5 * dxdt_sc(x_sc_arr[i - 1], v_sc_arr[i - 1])
				- 0.5 * dxdt_sc(x_sc_arr[i - 2], v_sc_arr[i - 2]));
		v_sc_arr[i] = v_sc_arr[i - 1]
			+ t_step_sc * (1.5 * dvdt_sc(x_sc_arr[i - 1], v_sc_arr[i - 1])
				- 0.5 * dvdt_sc(x_sc_arr[i - 2], v_sc_arr[i - 2]));
	}
	// scaled
	ofstream f_out(out_filename + "_sc.csv");
	f_out << "t,x,v" << endl;
	for (int i = 0; i < n; i++)
	{
		f_out << t_step_sc * i << ","
			<< x_sc_arr[i] << ","
			<< v_sc_arr[i] << endl;
	}
	f_out.close();

	// original
	f_out.open(out_filename + "_or.csv");
	f_out << "t,x,v" << endl;
	for (int i = 0; i < n; i++)
	{
		f_out << t_step * i << ","
			<< x_or(x_sc_arr[i]) << ","
			<< v_or(v_sc_arr[i]) << endl;
	}
	f_out.close();
	free(x_sc_arr);
	free(v_sc_arr);
}