#pragma once
#include "task2.h"
#include "my_data.h"
#include "matrix_funcs.h"
#include <vector>
#include <iostream>
#include <fstream>

#define ITER_LIM 100

using namespace std;

void t2(double x_step, double y_step, double& err, int& iter_num, int mode, std::string out_filename) {
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

	double w = 1.2;

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
				+f_func_1(x_step * i, y_step * j);
		}
	}

	for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++) {
			u_next[to_idx(i, j, N)] = -f_ext[to_idx(i, j, N)] / alpha;
		}
	}

	/////////////////////////////

	while (true) {
		u_prev = u_next;

		u_next[to_idx(1, 1, N)] =
			(1 - w) * u_prev[to_idx(1, 1, N)]
			- (f_ext[to_idx(1, 1, N)]
			+ beta * u_prev[to_idx(2, 1, N)]
			+ gamma * u_prev[to_idx(1, 2, N)]) * w / alpha;

		for (int i = 2; i < (N - 1); i++) {
			u_next[to_idx(i, 1, N)] =
				(1 - w) * u_prev[to_idx(i, 1, N)]
				- (f_ext[to_idx(i, 1, N)]
				+ beta * u_next[to_idx(i - 1, 1, N)]

				+ beta * u_prev[to_idx(i + 1, 1, N)]
				+ gamma * u_prev[to_idx(i, 2, N)]) * w / alpha;
		}

		u_next[to_idx((N - 1), 1, N)] =
			(1 - w) * u_prev[to_idx((N - 1), 1, N)]
			- (f_ext[to_idx(N - 1, 1, N)]
			+ beta * u_next[to_idx(N - 2, 1, N)]

			+ gamma * u_prev[to_idx(N - 1, 2, N)]) * w / alpha;


		for (int j = 2; j < (N - 1); j++) {
			u_next[to_idx(1, j, N)] =
				(1 - w) * u_prev[to_idx(1, j, N)]
				- (f_ext[to_idx(1, j, N)]
				+ gamma * u_next[to_idx(1, j - 1, N)]

				+ beta * u_prev[to_idx(2, j, N)]
				+ gamma * u_prev[to_idx(1, j + 1, N)]) * w / alpha;

			for (int i = 2; i < N - 1; i++) {
				u_next[to_idx(i, j, N)] =
					(1 - w) * u_prev[to_idx(i, j, N)]
					- (f_ext[to_idx(i, j, N)]
					+ gamma * u_next[to_idx(i, j - 1, N)]
					+ beta * u_next[to_idx(i - 1, j, N)]

					+ beta * u_prev[to_idx(i + 1, j, N)]
					+ gamma * u_prev[to_idx(i, j + 1, N)]) * w / alpha;
			}


			u_next[to_idx(N - 1, j, N)] =
				(1 - w) * u_prev[to_idx(N - 1, j, N)]
				- (f_ext[to_idx((N - 1), j, N)]
				+ gamma * u_next[to_idx(N - 1, j - 1, N)]
				+ beta * u_next[to_idx(N - 2, j, N)]

				+ gamma * u_prev[to_idx(N - 1, j + 1, N)]) * w / alpha;
		}

		u_next[to_idx(1, N - 1, N)] =
			(1 - w) * u_prev[to_idx(1, N - 1, N)]
			- (f_ext[to_idx(1, N - 1, N)]
			+ gamma * u_next[to_idx(1, N - 2, N)]

			+ beta * u_prev[to_idx(2, N - 1, N)]) * w / alpha;

		for (int i = 2; i < (N - 1); i++) {
			u_next[to_idx(i, N - 1, N)] =
				(1 - w) * u_prev[to_idx(i, N - 1, N)]
				- (f_ext[to_idx(i, N - 1, N)]
				+ gamma * u_next[to_idx(i, N - 2, N)]
				+ beta * u_next[to_idx(i - 1, N - 1, N)]

				+ beta * u_prev[to_idx(i + 1, N - 1, N)]) * w / alpha;
		}

		u_next[to_idx(N - 1, N - 1, N)] =
			(1 - w) * u_prev[to_idx(N - 1, N - 1, N)]
			- (f_ext[to_idx(N - 1, N - 1, N)]
			+ gamma * u_next[to_idx(N - 1, N - 2, N)]
			+ beta * u_next[to_idx(N - 2, N - 1, N)]) * w / alpha;

		/////////////////////////

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

	/////////////////////////////

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

void t2_mx(double x_step, double y_step, double& err, int& iter_num, int mode, std::string out_filename) {
	int N = LX / x_step;
	//int M = LY / y_step;   == N
	int iter_cnt = 0;
	double curr_err = 0;

	//double min_u_true = abs(u_func(x_step * 2, y_step * 2));
	vector<double> u_true((N - 1) * (N - 1));
	for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++) {
			u_true[to_idx(i, j, N)] = u_func(x_step * i, y_step * j);
			//min_u_true = abs(u_true[to_idx(i, j, N)]) min(min_u_true)
		}
	}

	double alpha = -2 * (1 / (x_step * x_step) + 1 / (y_step * y_step));
	double beta = 1 / (x_step * x_step);
	double gamma = 1 / (y_step * y_step);

	vector<vector<double>> A((N - 1) * (N - 1), vector<double>((N - 1) * (N - 1)));
	vector<vector<double>> B((N - 1) * (N - 1), vector<double>(1));
	
	for (int idx = 0; idx < (N - 1) * (N - 1); idx++) {
		A[idx][idx] = alpha;
		int i, j; to_ij(idx, N, i, j);
		B[idx][0] = -f_func_1(x_step * i, y_step * j);
	}
	for (int j = 1; j < N; j++) {
		for (int i = 1; i < N - 1; i++) {
			int idx = to_idx(i, j, N);
			A[idx][idx + 1] = beta;
			A[idx + 1][idx] = beta;
		}
	}
	for (int idx = 0; idx < (N - 1) * (N - 2); idx++) {
		A[idx][idx + (N - 1)] = gamma;
		A[idx + (N - 1)][idx] = gamma;
	}

	// i = 1 // i = N - 2 // j = 1 // j = N - 2
	for (int j = 1; j < N; j++) {
		int idx = to_idx(1, j, N);
		B[idx][0] += -phi_1(y_step * j) / (x_step * x_step);

		idx = to_idx(N - 1, j, N);
		B[idx][0] += -phi_2(y_step * j) / (x_step * x_step);

		idx = to_idx(j, 1, N);
		B[idx][0] += -phi_3(x_step * j) / (y_step * y_step);

		idx = to_idx(j, N - 1, N);
		B[idx][0] += -phi_4(x_step * j) / (y_step * y_step);
	}

	vector <vector <double>> u_sor((N - 1) * (N - 1), vector<double>(1));
	u_sor = B;
	u_sor = solve_sym_sor(A, B, err, 1.2, iter_num, mode);

	/*iter_num = 0;
	for (int i = 0; i < 5; i++) {
		u_sor = matrix_sub(u_sor, matrix_sub(matrix_prod(A, u_sor), B));
		iter_num++;
	}*/


	// save res
	ofstream f_out(out_filename + "_pars");
	f_out << N << endl << iter_num << endl << err << endl;
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
			f_out << u_sor[to_idx(i, j, N)][0] << " ";
		}
		f_out << endl;
	}
	f_out.close();
}