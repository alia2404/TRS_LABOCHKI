#include <iostream>
#include <fstream>
#include "task5.h"

using namespace std;

void t5(int n, double x_start, double x_fin, string out_filename, double& x_step) {
	// n - number of parts; number of nodes = n + 1
	x_step = (x_fin - x_start) / n;
	shotgun_method(n, x_start, x_fin, out_filename);
}

void shotgun_method(int n, double x_start, double x_fin, string out_filename) {
	double teta1 = 10, teta2 = 50;
	double psi1 = rk4_search(n, x_start, x_fin, teta1, "");
	double psi2 = rk4_search(n, x_start, x_fin, teta2, "");
	printf("psi1 = %lf,\t psi2 = %lf\n", psi1, psi2);
	double teta0 = (psi2 * teta1 - psi1 * teta2) / (psi2 - psi1);
	double psi0 = rk4_search(n, x_start, x_fin, teta0, out_filename);
	printf("psi0 = %lf\n", psi0);
}

double dudt(double x, double u, double v) { return v; }
double dvdt(double x, double u, double v) { return u + 2 - x*x; }

double rk4_search(int n, double x_start, double x_fin, double u_start, string out_filename) {
	double x_step = (x_fin - x_start) / n;
	double* u; u = (double*)malloc((n + 1) * sizeof(double));
	double* v; v = (double*)malloc((n + 1) * sizeof(double));
	u[0] = u_start;
	v[0] = u_start - 1;
	for (int i = 0; i < n; i++) {
		double x_i = x_start + x_step * i;
		double k1_u = dudt(x_i, u[i], v[i]);
		double k1_v = dvdt(x_i, u[i], v[i]);
		double k2_u = dudt(x_i + 0.5 * x_step, u[i] + 0.5 * x_step * k1_u, v[i] + 0.5 * x_step * k1_v);
		double k2_v = dvdt(x_i + 0.5 * x_step, u[i] + 0.5 * x_step * k1_u, v[i] + 0.5 * x_step * k1_v);
		double k3_u = dudt(x_i + 0.5 * x_step, u[i] + 0.5 * x_step * k2_u, v[i] + 0.5 * x_step * k2_v);
		double k3_v = dvdt(x_i + 0.5 * x_step, u[i] + 0.5 * x_step * k2_u, v[i] + 0.5 * x_step * k2_v);
		double k4_u = dudt(x_i + x_step, u[i] + x_step * k3_u, v[i] + x_step * k3_v);
		double k4_v = dvdt(x_i + x_step, u[i] + x_step * k3_u, v[i] + x_step * k3_v);
		u[i + 1] = u[i] + x_step * (k1_u + 2 * k2_u + 2 * k3_u + k4_u) / 6;
		v[i + 1] = v[i] + x_step * (k1_v + 2 * k2_v + 2 * k3_v + k4_v) / 6;
		if (i % 1000 == 0)
			printf("%lf %lf, \t%lf %lf, \t%lf %lf, \t%lf %lf\n", k1_u, k1_v, k2_u, k2_v, k3_u, k3_v, k4_u, k4_v);
	}
	if (out_filename != "") {
		ofstream f_out(out_filename);
		f_out << "x,u" << endl;
		for (int i = 0; i <= n; i++) {
			f_out << x_start + x_step * i << "," << u[i] << endl;
		}
		f_out.close();
	}
	printf("u[n-1] = %lf,\t u[n] = %lf\n", u[n - 1], u[n]);
	printf("v[n-1] = %lf,\t v[n] = %lf\n", v[n - 1], v[n]);
	double v_n = v[n];
	free(u);
	free(v);
	return psi(v_n);
}

double psi(double v) {
	return v - 1;
}