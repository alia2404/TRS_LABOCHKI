#pragma once

#include "task4.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <chrono>

using namespace std;

void t4(double x_step, double t_step, double t_fin, double sigma, std::string out_filename) {
	if (sigma < 0.00001) {
		explicit_t4(x_step, t_step, t_fin, out_filename);
	}
	else
	{	
		implicit_t4(x_step, t_step, t_fin, sigma, out_filename);
	}
}

void explicit_t4(double x_step, double t_step, double t_fin, std::string out_filename) {
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
	for (int j = 0; j < L; j++) {
		// calc inner nodes
		for (int i = 1; i < N; i++) {
			double t_j = t_step * j;
			double x_i = x_step * i;
			// cout << i + (j + 1) * (L + 1) << endl;
			double k_w = k((u[j][i] + u[j][i - 1]) / 2);
			double k_e = k((u[j][i] + u[j][i + 1]) / 2);
			u[j + 1][i] =
				(1 - (k_e + k_w) * t_step / (x_step * x_step)) * u[j][i]
				+ t_step * (u[j][i + 1] * k_w + u[j][i - 1] * k_e) / (x_step * x_step)
				+ Func(u[j][i]) * func(x_i, t_j) * t_step;
		}
		// calc border nodes using border conditions
		// u_0^j - ux_0^j = psi1(t_j)
		// ux_N^j = psi2(t_j)
		u[j + 1][0] = (psi0(t_step * (j + 1)) + u[j + 1][1] / x_step) * x_step / (1 + x_step);

		u[j + 1][N] = psi1(t_step * (j + 1)) * x_step + u[j + 1][N - 1];
	}

	auto end_time = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> elapsed_time = end_time - start_time;
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

	ofstream f_out(out_filename + "_pars");
	f_out << N << " " << L << " " << x_step << " " << t_step << " " << 0. << endl;
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

void implicit_t4(double x_step, double t_step, double t_fin, double sigma, std::string out_filename) {
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
	for (int j = 0; j < L; j++) {
		// 3diag gauss
		// fill matrix
		vector<double> a_coef(N + 1);
		vector<double> b_coef(N + 1);
		vector<double> c_coef(N + 1);
		vector<double> d_coef(N + 1);
		a_coef[0] = 0;
		b_coef[0] = 1 + 1 * x_step;
		c_coef[0] = 1;
		d_coef[0] = -x_step * psi0(t_step * (j + 1));

		a_coef[N] = 1;
		b_coef[N] = 1;
		c_coef[N] = 0;
		d_coef[N] = -x_step * psi1(t_step * (j + 1));

		for (int i = 1; i < N; i++) {
			double t_j = t_step * j;
			double x_i = x_step * i;
			//double k_w_1 = k((u[j + 1][i] + u[j + 1][i + 1]) / 2);
			//double k_e_1 = k((u[j + 1][i] + u[j + 1][i - 1]) / 2);
			double k_w = k((u[j][i] + u[j][i - 1]) / 2);
			double k_e = k((u[j][i] + u[j][i + 1]) / 2);

			a_coef[i] = sigma * k_w * t_step / (x_step * x_step);
			b_coef[i] = sigma*(k_e + k_w) * t_step / (x_step * x_step) + 1;
			c_coef[i] = sigma * k_e * t_step / (x_step * x_step);
			d_coef[i] = -t_step * Func(u[j][i]) * func(x_i, t_j)
				- (1 - sigma) * (k_e * (u[j][i + 1] - u[j][i]) - k_w * (u[j][i] - u[j][i - 1])) * t_step / (x_step * x_step)
				- u[j][i];
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
	f_out << N << " " << L << " " << x_step << " " << t_step << " " << sigma << endl;
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