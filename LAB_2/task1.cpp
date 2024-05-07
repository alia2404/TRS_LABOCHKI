#pragma once

#include "task1.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <chrono>

using namespace std;

void t1(double x_step, double t_step, double t_fin, std::string out_filename) {
	int N = 1 / x_step;
	int L = t_fin / t_step;
	vector<vector<double>> u = vector<vector<double>>(L+1, vector<double>(N + 1));
	printf("Vector len: %d x %d\n", u.size(), u[0].size());

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
			// t_j = t_step * j
			// x_i = x_step * i
			// cout << i + (j + 1) * (L + 1) << endl;
			u[j + 1][i] =
				(1 - 2 * t_step / (x_step * x_step)) * u[j][i]
				+ t_step / (x_step * x_step) * (u[j][i - 1] + u[j][i + 1])
				+ t_step * func(x_step * i, t_step * j);
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