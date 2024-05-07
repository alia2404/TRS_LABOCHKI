#include <iostream>
#include <fstream>
#include "task3.h"

using namespace std;

void t3(double t_step, double t_fin, string out_filename) {
	int n = t_fin / t_step + 1;
	double t_step_sc = t_sc(t_step);
	double t_fin_sc = t_sc(t_fin);
	double* x_sc_arr; x_sc_arr = (double*) malloc(n * sizeof(double));
	double* v_sc_arr; v_sc_arr = (double*) malloc(n * sizeof(double));
	x_sc_arr[0] = x_sc(x0);
	v_sc_arr[0] = v_sc(v0);
	//x_sc_arr[1] = x_sc_arr[0] + t_step_sc * dxdt_sc(x_sc_arr[0], v_sc_arr[0]);
	//v_sc_arr[1] = v_sc_arr[0] + t_step_sc * dvdt_sc(x_sc_arr[0], v_sc_arr[0]);
	for (int i = 1; i < n; i++) {
		double k1_x = dxdt_sc(x_sc_arr[i - 1], v_sc_arr[i - 1]);
		double k1_v = dvdt_sc(x_sc_arr[i - 1], v_sc_arr[i - 1]);
		double k2_x = dxdt_sc(x_sc_arr[i - 1] + 0.5*t_step_sc*k1_x,
			v_sc_arr[i - 1] + 0.5 * t_step_sc * k1_v);
		double k2_v = dvdt_sc(x_sc_arr[i - 1] + 0.5 * t_step_sc * k1_x,
			v_sc_arr[i - 1] + 0.5 * t_step_sc * k1_v);
		double k3_x = dxdt_sc(x_sc_arr[i - 1] + 0.5 * t_step_sc * k2_x,
			v_sc_arr[i - 1] + 0.5 * t_step_sc * k2_v);
		double k3_v = dvdt_sc(x_sc_arr[i - 1] + 0.5 * t_step_sc * k2_x,
			v_sc_arr[i - 1] + 0.5 * t_step_sc * k2_v);
		double k4_x = dxdt_sc(x_sc_arr[i - 1] + t_step_sc * k3_x,
			v_sc_arr[i - 1] + t_step_sc * k3_v);
		double k4_v = dvdt_sc(x_sc_arr[i - 1] + t_step_sc * k3_x,
			v_sc_arr[i - 1] + t_step_sc * k3_v);
		x_sc_arr[i] = x_sc_arr[i - 1]
			+ t_step_sc * (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6;
		v_sc_arr[i] = v_sc_arr[i - 1]
			+ t_step_sc * (k1_v + 2 * k2_v + 2 * k3_v + k4_v) / 6;
		//if (i % 1000 == 0)
		//	printf("%lf %lf, \t%lf %lf, \t%lf %lf, \t%lf %lf\n", k1_x, k1_v, k2_x, k2_v, k3_x, k3_v, k4_x, k4_v);
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