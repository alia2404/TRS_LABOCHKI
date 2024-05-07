#pragma once

#include "task5.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <chrono>

using namespace std; 

void t5(double x_step, double t_step, double t_fin, double delta, std::string out_filename) {
	int N = 1 / x_step;
	int L = t_fin / t_step;
	vector<int> num_of_iters(L + 1, 0);
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

	vector<double> u_prev = vector<double>(N + 1);
	vector<double> u_next = vector<double>(N + 1);

	vector<double> k_w_prev = vector<double>(N - 1, 0);
	vector<double> k_w_next = vector<double>(N - 1, 0);

	vector<double> k_e_prev = vector<double>(N - 1, 0);
	vector<double> k_e_next = vector<double>(N - 1, 0);

	double delta_k_e = 2 * delta;
	double delta_k_w = 2 * delta;

	for (int j = 0; j < L; j++) {
		for (int i = 0; i < N + 1; i++) {
			u_next[i] = u[j][i];
		}

		int iter_count = 0;
		do {
			u_prev = u_next;
			k_w_prev = k_w_next;
			k_e_prev = k_e_next;

			d_coef[0] = -x_step * psi0(t_step * (j + 1));
			d_coef[N] = -x_step * psi1(t_step * (j + 1));

			for (int i = 1; i < N; i++) {
				double t_j = t_step * j;
				double x_i = x_step * i;
				double k_w = k((u_prev[i] + u_prev[i - 1]) / 2);
				double k_e = k((u_prev[i] + u_prev[i + 1]) / 2);

				k_w_next[i - 1] = k_w;
				k_e_next[i - 1] = k_e;

				a_coef[i] = k_w * t_step / (x_step * x_step);
				b_coef[i] = (k_e + k_w) * t_step / (x_step * x_step) + 1;
				c_coef[i] = k_e * t_step / (x_step * x_step);
				d_coef[i] = -t_step * Func(u[j][i]) * func(x_i, t_j) - u[j][i];
			}

			diag3_method(a_coef, b_coef, c_coef, d_coef, N + 1, u_next);

			iter_count++;
			if (iter_count < 2) continue;

			delta_k_e = calc_delta(k_e_prev, k_e_next);
			delta_k_w = calc_delta(k_e_prev, k_e_next);

			if (delta > delta_k_e && delta > delta_k_w) break;

		} while (true);

		num_of_iters[j + 1] = iter_count;
		u[j + 1] = u_next;

	}
	auto end_time = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> elapsed_time = end_time - start_time;
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

	ofstream f_out(out_filename + "_pars");
	f_out << N << " " << L << " " << x_step << " " << t_step << " " << delta << endl;
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
	f_out.open(out_filename + "_iters.csv");
	f_out << "t,iters" << endl;
	for (int j = 0; j < L + 1; j++) {
		f_out << t_step * j << "," << num_of_iters[j] << endl;
	}
	f_out.close();
}

double calc_delta(vector<double> v1, vector<double> v2) {
	double delta = abs(v1[0] - v2[0]);
	for (int i = 1; i < v2.size(); i++) {
		delta = max(delta, abs(v1[i] - v2[i]));
	}
	return delta;
}


void t5_u(double x_step, double t_step, double t_fin, double delta, std::string out_filename) {
	int N = 1 / x_step;
	int L = t_fin / t_step;
	vector<int> num_of_iters(L + 1, 0);
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

	for (int j = 0; j < L; j++) {
		vector<double> u_prev = vector<double>(N + 1);
		vector<double> u_next = vector<double>(N + 1);
		for (int i = 0; i < N + 1; i++) {
			u_next[i] = u[j][i];
		}

		int iter_count = 0;
		do {
			u_prev = u_next;
			d_coef[0] = -x_step * psi0(t_step * (j + 1));
			d_coef[N] = -x_step * psi1(t_step * (j + 1));

			for (int i = 1; i < N; i++) {
				double t_j = t_step * j;
				double x_i = x_step * i;
				double k_w = k((u_prev[i] + u_prev[i - 1]) / 2);
				double k_e = k((u_prev[i] + u_prev[i + 1]) / 2);

				a_coef[i] = k_w * t_step / (x_step * x_step);
				b_coef[i] = (k_e + k_w) * t_step / (x_step * x_step) + 1;
				c_coef[i] = k_e * t_step / (x_step * x_step);
				d_coef[i] = -t_step * Func(u[j][i]) * func(x_i, t_j) - u[j][i];
			}

			diag3_method(a_coef, b_coef, c_coef, d_coef, N + 1, u_next);

			iter_count++;
		} while (delta < calc_delta(u_prev, u_next));
		num_of_iters[j + 1] = iter_count;
		u[j + 1] = u_next;

	}


	auto end_time = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> elapsed_time = end_time - start_time;
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

	ofstream f_out(out_filename + "_pars");
	f_out << N << " " << L << " " << x_step << " " << t_step << " " << delta << endl;
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
	f_out.open(out_filename + "_iters.csv");
	f_out << "t,iters" << endl;
	for (int j = 0; j < L + 1; j++) {
			f_out << t_step * j << "," << num_of_iters[j] << endl;
	}
	f_out.close();
}


double calc_delta_rel(vector<double> v1, vector<double> v2) {
	double delta = abs((v1[0] - v2[0]) / v2[0]);
	for (int i = 1; i < v2.size(); i++) {
		delta = max(delta, abs((v1[i] - v2[i]) / v2[i]));
	}
	return delta;
}