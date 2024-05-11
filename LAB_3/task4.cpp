#pragma once
#include "task4.h"
#include "my_data.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "eigen-3.4.0\Eigen\Eigen.h"

#define ITER_LIM 100

using namespace std;

void t4(double x_step, double y_step, double& err, int& iter_num, int mode, std::string out_filename) {
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

	// size of the image
	// number of unknowns (=number of pixels)
	
	std::vector<Trip> coefficients;            // list of non-zeros coefficients
	Eigen::VectorXd B((N - 1) * (N - 1));                  // the right hand side-vector resulting from the constraints
	//buildProblem(coefficients, b, n);
	for (int idx = 0; idx < (N - 1) * (N - 1); idx++) {
		int i, j; to_ij(idx, N, i, j);
		B[idx] = -f_func_2(x_step * i, y_step * j);
	}

	for (int j = 1; j < N; j++) {
		for (int i = 1; i < N; i++) {
			int idx = to_idx(i, j, N);
			double alpha = c_func(x_step * i, y_step * j)
				- a_func(x_step * (i + 0.5), y_step * j) / (x_step * x_step)
				- a_func(x_step * (i - 0.5), y_step * j) / (x_step * x_step)
				- b_func(x_step * i, y_step * (j + 0.5)) / (y_step * y_step)
				- b_func(x_step * i, y_step * (j - 0.5)) / (y_step * y_step);
			coefficients.push_back(Trip(idx, idx, alpha));
		}
		for (int i = 2; i < N; i++) {
			int idx = to_idx(i, j, N);
			double beta = a_func(x_step * (i - 0.5), y_step * j) / (x_step * x_step);
			coefficients.push_back(Trip(idx, idx - 1, beta));
		}
		for (int i = 1; i < N - 1; i++) {
			int idx = to_idx(i, j, N);
			double beta = a_func(x_step * (i + 0.5), y_step * j) / (x_step * x_step);
			coefficients.push_back(Trip(idx, idx + 1, beta));
		}
	}
	for (int j = 2; j < N; j++) {
		for (int i = 1; i < N; i++) {
			int idx = to_idx(i, j, N);
			double gamma = b_func(x_step * i, y_step * (j - 0.5)) / (y_step * y_step);
			coefficients.push_back(Trip(idx, idx - (N - 1), gamma));
		}
	}
	for (int j = 1; j < N - 1; j++) {
		for (int i = 1; i < N; i++) {
			int idx = to_idx(i, j, N);
			double gamma = b_func(x_step * i, y_step * (j + 0.5)) / (y_step * y_step);
			coefficients.push_back(Trip(idx, idx + (N - 1), gamma));
		}
	}

	// i = 1 // i = N - 2 // j = 1 // j = N - 2
	for (int j = 1; j < N; j++) {
		int idx = to_idx(1, j, N);
		B[idx] += -phi_1(y_step * j) / (x_step * x_step);

		idx = to_idx(N - 1, j, N);
		B[idx] += -phi_2(y_step * j) / (x_step * x_step);

		idx = to_idx(j, 1, N);
		B[idx] += -phi_3(x_step * j) / (y_step * y_step);

		idx = to_idx(j, N - 1, N);
		B[idx] += -phi_4(x_step * j) / (y_step * y_step);
	}


	SpMat A((N - 1) * (N - 1), (N - 1) * (N - 1));
	A.setFromTriplets(coefficients.begin(), coefficients.end());
	/*
	for (int i = 0; i < (N-2)*(N-2); i++) {
		for (int j = 0; j < (N - 2) * (N - 2); j++) {
			cout << A.coeff(i, j) << " ";
		}
		cout << endl;
	}
	for (int i = 0; i < (N - 2) * (N - 2); i++) {
		cout << B.coeff(i) << " ";
		cout << endl;
	}

	cout << "----------" << endl;
	
	for (int j = 1; j < N - 1; j++) {
		for (int i = 1; i < N - 1; i++) {
			cout << -f_func_1(x_step * i, y_step * j) << endl;
		}
	}
	*/

	// Solving
	Eigen::SimplicialCholesky<SpMat> chol(A);  // performs a Cholesky factorization of A
	Eigen::VectorXd u_eig = chol.solve(B);         // use the factorization to solve for the given right hand side

	// Export the result to a file:
	//saveAsBitmap(x, n, argv[1]);

	///////////////////////////////////////////


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
			f_out << u_eig[to_idx(i, j, N)] << " ";
		}
		f_out << endl;
	}
	f_out.close();
}