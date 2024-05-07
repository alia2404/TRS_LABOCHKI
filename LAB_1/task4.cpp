#include <iostream>
#include <fstream>
#include "task4.h"

using namespace std;

void t4(int n, double x_start, double x_fin, string out_filename, double& x_step) {
	// n - number of parts; number of nodes = n + 1
	x_step = (x_fin - x_start) / n;
	double h = x_step;
	double* a = (double*)malloc((n + 1) * sizeof(double));
	double* b = (double*)malloc((n + 1) * sizeof(double));
	double* c = (double*)malloc((n + 1) * sizeof(double));
	double* d = (double*)malloc((n + 1) * sizeof(double));
	a[0] = 0; b[0] = -(1 + 1 / h + h / 2); c[0] = -1 / h;
	d[0] = 1 + (x_start * x_start - 2) * h / 2;
	a[n] = -1 / h; b[n] = -(1 / h + h / 2); c[n] = 0;
	d[n] = 1 + (x_fin * x_fin - 2) * h / 2;
	for (int i = 1; i < n; i++) {
		a[i] = 1;
		b[i] = h * h + 2;
		c[i] = 1;
		double x_i = x_start + i * h;
		d[i] = (2 - x_i * x_i) * (h * h);
	}

	double* u = (double*)malloc((n + 1) * sizeof(double));
	diag3_method(a, b, c, d, n + 1, u);
	free(a); free(b); free(c); free(d);
	ofstream f_out(out_filename);
	f_out << "x,u" << endl;
	for (int i = 0; i <= n; i++) {
		f_out << x_start + x_step * i << "," << u[i] << endl;
	}
	f_out.close();
	free(u);
}

void diag3_method(double* a, double* b, double* c, double* d, int n, double* x) {
	/*
	* a_i*x_{i-1} - b_i*x_i + c_i*x_{i+1} = d_i
	* a_0 = c_{n-1} = 0
	* i = 0..(n-1)
	* n equations
	*/
	double* ksi = (double*)malloc((n + 1) * sizeof(double));
	double* teta = (double*)malloc((n + 1) * sizeof(double));
	ksi[0] = 0;
	teta[0] = 0;
	for (int i = 0; i < n; i++) {
		double denominator = b[i] - a[i] * ksi[i];
		ksi[i + 1] = c[i] / denominator;
		teta[i + 1] = (a[i]*teta[i] - d[i]) / denominator;
	}
	free(x);
	x = (double*)malloc(n * sizeof(double));
	x[n - 1] = teta[n];
	for (int i = n-2; i >= 0; i--) {
		x[i] = ksi[i + 1] * x[i + 1] + teta[i + 1];
	}
	free(ksi);
	free(teta);
}