#pragma once
#include "task1.h"
#include "my_data.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "eigen-3.4.0\Eigen\Core"


using namespace std;

void t1(double x_step, double y_step, double& err, int& iter_num, int mode, std::string out_filename) {
	int N = LX / x_step;
	//int M = LY / y_step; == N
	vector<double> u_prev((N - 1) * (N - 1), phi_1(0));
	vector<double> u_next((N - 1) * (N - 1), phi_1(0));
	vector<double> f_ext((N - 1) * (N - 1));
	int iter_cnt = 0;
	double curr_err = 0;

	vector<double> u_true((N - 1) * (N - 1));
	for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++) {
			u_true[to_idx(i, j, N)] = u_func(x_step * i, y_step * j);
		}
	}

	double alpha = -2 * (1 / (x_step * x_step) + 1 / (y_step * y_step));
	double beta = 1 / (x_step * x_step);
	double gamma = 1 / (y_step * y_step);

	f_ext[to_idx(1, 1, N)] =
		+f_func_1(x_step * 1, y_step * 1)
		+ phi_1(y_step * 1) / (x_step * x_step)
		+ phi_3(x_step * 1) / (y_step * y_step);

	f_ext[to_idx(1, N - 1, N)] =
		+f_func_1(x_step * 1, y_step * (N - 1))
		+ phi_1(y_step * (N - 1)) / (x_step * x_step)
		+ phi_4(x_step * 1) / (y_step * y_step);

	f_ext[to_idx(N - 1, 1, N)] =
		+f_func_1(x_step * (N - 1), y_step * 1)
		+ phi_2(y_step * 1) / (x_step * x_step)
		+ phi_3(x_step * (N - 1)) / (y_step * y_step);

	f_ext[to_idx(N - 1, N - 1, N)] =
		+f_func_1(x_step * (N - 1), y_step * (N - 1))
		+ phi_2(y_step * (N - 1)) / (x_step * x_step)
		+ phi_4(x_step * (N - 1)) / (y_step * y_step);

	for (int i = 2; i < (N - 1); i++) {
		f_ext[to_idx(i, 1, N)] =
			+f_func_1(x_step * i, y_step * 1)
			+ phi_3(x_step * i) / (y_step * y_step);
		f_ext[to_idx(i, N - 1, N)] =
			+f_func_1(x_step * i, y_step * (N - 1))
			+ phi_4(x_step * i) / (y_step * y_step);
	}

	for (int j = 2; j < (N - 1); j++) {
		f_ext[to_idx(1, j, N)] =
			+f_func_1(x_step * 1, y_step * j)
			+ phi_1(y_step * j) / (x_step * x_step);
		f_ext[to_idx(N - 1, j, N)] =
			+f_func_1(x_step * (N - 1), y_step * j)
			+ phi_2(y_step * j) / (x_step * x_step);
	}

	for (int i = 2; i < N - 1; i++) {
		for (int j = 2; j < N - 1; j++) {
			f_ext[to_idx(i, j, N)] =
				+ f_func_1(x_step * i, y_step * j);
		}
	}

	for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++) {
			u_next[to_idx(i, j, N)] = -f_ext[to_idx(i, j, N)] / alpha;
		}
	}

	while (true) {
		u_prev = u_next;

		u_next[to_idx(1, 1, N)] =
			-(f_ext[to_idx(1, 1, N)]
			+ beta * u_prev[to_idx(2, 1, N)]
			+ gamma * u_prev[to_idx(1, 2, N)]) / alpha;

		u_next[to_idx(1, N - 1, N)] =
			-(f_ext[to_idx(1, N - 1, N)]
			+ gamma * u_prev[to_idx(1, N - 2, N)]
			+ beta * u_prev[to_idx(2, N - 1, N)]) / alpha;
			
		u_next[to_idx((N - 1), 1, N)] =
			-(f_ext[to_idx(N - 1, 1, N)]
			+ beta * u_prev[to_idx(N - 2, 1, N)]
			+ gamma * u_prev[to_idx(N - 1, 2, N)]) / alpha;

		u_next[to_idx(N - 1, N - 1, N)] =
			-(f_ext[to_idx(N - 1, N - 1, N)]
			+ gamma * u_prev[to_idx(N - 1, N - 2, N)]
			+ beta * u_prev[to_idx(N - 2, N - 1, N)]) / alpha;

		for (int i = 2; i < (N - 1); i++) {
			u_next[to_idx(i, 1, N)] =
				-(f_ext[to_idx(i, 1, N)]
				+ beta * u_prev[to_idx(i - 1, 1, N)]
				+ beta * u_prev[to_idx(i + 1, 1, N)]
				+ gamma * u_prev[to_idx(i, 2, N)]) / alpha;

			u_next[to_idx(i, N - 1, N)] =
				-(f_ext[to_idx(i, N - 1, N)]
				+ gamma * u_prev[to_idx(i, N - 2, N)]
				+ beta * u_prev[to_idx(i - 1, N - 1, N)]
				+ beta * u_prev[to_idx(i + 1, N - 1, N)]) / alpha;

		}

		for (int j = 2; j < (N - 1); j++) {
			u_next[to_idx(1, j, N)] =
				-(f_ext[to_idx(1, j, N)]
				+ gamma * u_prev[to_idx(1, j - 1, N)]
				+ beta * u_prev[to_idx(2, j, N)]
				+ gamma * u_prev[to_idx(1, j + 1, N)]) / alpha;

			u_next[to_idx(N - 1, j, N)] =
				-(f_ext[to_idx((N - 1), j, N)]
				+ gamma * u_prev[to_idx(N - 1, j - 1, N)]
				+ beta * u_prev[to_idx(N - 2, j, N)]
				+ gamma * u_prev[to_idx(N - 1, j + 1, N)]) / alpha;
				
		}

		for (int i = 2; i < N - 1; i++) {
			for (int j = 2; j < N - 1; j++) {
				u_next[to_idx(i, j, N)] =
					-(f_ext[to_idx(i, j, N)]
					+ gamma * u_prev[to_idx(i, j - 1, N)]
					+ beta * u_prev[to_idx(i - 1, j, N)]
					+ beta * u_prev[to_idx(i + 1, j, N)]
					+ gamma * u_prev[to_idx(i, j + 1, N)]) / alpha;
			}
		}

		iter_cnt++;
		if (mode == ITER_CONST) {
			if (iter_cnt >= iter_num) {
				err = calc_delta(u_next, u_prev);
				break;
			}
		} 
		else /*if (mode == ERR_CONST)*/ {
			curr_err = calc_delta(u_next, u_prev);
			if (curr_err <= err) {
				iter_num = iter_cnt;
				break;
			}
		}
		if (iter_cnt % 100 == 0) cout << "Iter number " << iter_cnt << "; max_diff = " << curr_err << endl;
	}

	// save res
	ofstream f_out(out_filename + "_pars");
	f_out << N << endl << iter_cnt << endl << err << endl;
	f_out.close();

	f_out.open(out_filename + "_true_vals");
	for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++) {
			f_out << u_true[to_idx(i, j, N)] << " ";
		}
		f_out << endl;
	}
	f_out.close();

	f_out.open(out_filename + "_vals");
	for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++) {
			f_out << u_next[to_idx(i, j, N)] << " ";
		}
		f_out << endl;
	}
	f_out.close();
}