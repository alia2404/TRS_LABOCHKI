#pragma once

#include "task3.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <chrono>

using namespace std;

void t3(double x_step, double t_step, double t_fin, std::string out_filename) {
	int N = 1 / x_step;
	int L = t_fin / t_step;
	vector<vector<double>> u = vector<vector<double>>(L + 1, vector<double>(N + 1));
	printf("Vector len: %d\n", u.size());

	auto start_time = std::chrono::high_resolution_clock::now();
	// fill zero layer
	for (int i = 0; i < N + 1; i++) {
		// t = 0; x = i * x_step
		u[0][i] = phi(i * x_step);
	}
	// calc layers
	vector<double> a_coef(N + 1);
	vector<double> b_coef(N + 1);
	vector<double> c_coef(N + 1);
	vector<double> d_coef(N + 1);
	a_coef[0] = 0;
	b_coef[0] = 1 + 1 * x_step;
	c_coef[0] = 1;

	a_coef[N] = 1;
	b_coef[N] = 1;
	c_coef[N] = 0;

	for (int i = 1; i < N; i++) {
		a_coef[i] = 1;       
		b_coef[i] = 2 + 2 * x_step * x_step / t_step;
		c_coef[i] = 1;
	}
	for (int j = 0; j < L; j++) {
		// 3diag gauss
		// fill matrix
		d_coef[0] = -x_step * psi0(t_step * (j + 1));
		d_coef[N] = -x_step * psi1(t_step * (j + 1));

		for (int i = 1; i < N; i++) {
			d_coef[i] = -(u[j][i - 1] + u[j][i + 1])
				- 2 * x_step * x_step * func(i * x_step, j * t_step)
				+ 2 * u[j][i] * (1 - x_step * x_step / t_step);
			//-u[j][i] * x_step * x_step / t_step - x_step * x_step * func(x_i, j * t_step);
		}

		vector<double> u_next(N + 1);
		diag3_method(a_coef, b_coef, c_coef, d_coef, N + 1, u_next);

		u[j + 1] = u_next;
		//printf("%lf%% is processed\n", (j + 1) * 100. / L);
	}

	auto end_time = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> elapsed_time = end_time - start_time;
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

	ofstream f_out(out_filename + "_pars");
	f_out << N << " " << L << " " << x_step << " " << t_step << endl;
	f_out.close();

	f_out.open(out_filename + "_time");
	cout << elapsed_time.count() << " ms" << endl;
	f_out << elapsed_time.count() << endl;
	f_out.close();

	f_out.open(out_filename + "_vals.csv");
	f_out << "t,x,u" << endl;
	for (int j = 0; j < L + 1; j++) {
		for (int i = 0; i < N + 1; i++) {
			f_out << t_step * j << ","
				<< x_step * i << ","
				<< std::fixed << u[j][i] << endl;
		}
	}
	f_out.close();
}