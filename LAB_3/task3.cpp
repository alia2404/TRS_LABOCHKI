#pragma once
#include "task3.h"
#include "my_data.h"
#include "matrix_funcs.h"
#include <vector>
#include <iostream>
#include <fstream>

#define ITER_LIM 100

using namespace std;

void t3(double x_step, double y_step, double& err, int& iter_num, int mode, std::string out_filename) {
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
	
	//std::vector<Trip> coefficients;            // list of non-zeros coefficients
	//Eigen::VectorXd B((N - 2) * (N - 2));                  // the right hand side-vector
	vector<vector<double>> A((N - 1) * (N - 1), vector<double>((N - 1) * (N - 1), 0));                  // the right hand side-vector
	vector<double> B((N - 1) * (N - 1), 0);                  // the right hand side-vector

	for (int idx = 0; idx < (N - 1) * (N - 1); idx++) {
		int i, j; to_ij(idx, N, i, j);
		//coefficients.push_back(Trip(idx, idx, alpha));
		A[idx][idx] = alpha;
		B[idx] = -f_func_1(x_step * i, y_step * j);
	}
	for (int j = 1; j < N; j++) {
		for (int i = 1; i < N - 1; i++) {
			int idx = to_idx(i, j, N);
			//coefficients.push_back(Trip(idx, idx + 1, beta));
			//coef2.push_back(Trip(idx, idx+1, subdiag_val));
			A[idx][idx + 1] = beta;

			//coefficients.push_back(Trip(idx + 1, idx, beta));
			//coef2.push_back(Trip(idx+1, idx, subdiag_val));
			A[idx + 1][idx] = beta;
		}
	}
	for (int idx = 0; idx < (N - 1) * (N - 2); idx++) {
		//coefficients.push_back(Trip(idx, idx + (N - 2), gamma));
		//coef2.push_back(Trip(idx, idx + (N - 2), fardiag_val));
		A[idx][idx + (N - 1)] = gamma;

		//coefficients.push_back(Trip(idx + (N - 2), idx, gamma));
		//coef2.push_back(Trip(idx + (N - 2), idx, fardiag_val));
		A[idx + (N - 1)][idx] = gamma;
	}
	

	// i = 1 // i = N - 1 // j = 1 // j = N - 1
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


	//SpMat A((N - 2) * (N - 2), (N - 2) * (N - 2));
	//A.setFromTriplets(coefficients.begin(), coefficients.end());

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
	vector<double> u_sq((N - 1) * (N - 1));
	LU_decay(A, B, u_sq);

	cout << "finish" << endl;

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
			f_out << u_sq[to_idx(i, j, N)] << " ";
		}
		f_out << endl;
	}
	f_out.close();
}


void t3_eigen(double x_step, double y_step, double& err, int& iter_num, int mode, std::string out_filename) {
	int N = LX / x_step + 1;
	//int M = LY / y_step + 1;   == N
	int iter_cnt = 0;
	double curr_err = 0;

	//double min_u_true = abs(u_func(x_step * 2, y_step * 2));
	vector<double> u_true((N - 2) * (N - 2));
	for (int i = 1; i < N - 1; i++) {
		for (int j = 1; j < N - 1; j++) {
			u_true[to_idx(i, j, N)] = u_func(x_step * i, y_step * j);
			//min_u_true = abs(u_true[to_idx(i, j, N)]) min(min_u_true)
		}
	}

	double alpha = -2 * (1 / (x_step * x_step) + 1 / (y_step * y_step));
	double beta = 1 / (x_step * x_step);
	double gamma = 1 / (y_step * y_step);

	std::vector<Trip> coefficients;            // list of non-zeros coefficients
	Eigen::VectorXd B((N - 2) * (N - 2));                  // the right hand side-vector

	std::vector<Trip> coef2;
	double diag_val = -0.1;
	double subdiag_val = 0;
	double fardiag_val = 0;


	for (int idx = 0; idx < (N - 2) * (N - 2); idx++) {
		int i, j; to_ij(idx, N, i, j);
		coefficients.push_back(Trip(idx, idx, alpha));
		coef2.push_back(Trip(idx, idx, diag_val));
		B[idx] = -f_func_1(x_step * i, y_step * j);
	}
	for (int j = 1; j < N - 1; j++) {
		for (int i = 1; i < N - 2; i++) {
			int idx = to_idx(i, j, N);
			coefficients.push_back(Trip(idx, idx + 1, beta));
			coef2.push_back(Trip(idx, idx + 1, subdiag_val));
			coefficients.push_back(Trip(idx + 1, idx, beta));
			coef2.push_back(Trip(idx + 1, idx, subdiag_val));
		}
	}
	for (int idx = 0; idx < (N - 2) * (N - 3); idx++) {
		coefficients.push_back(Trip(idx, idx + (N - 2), gamma));
		coef2.push_back(Trip(idx, idx + (N - 2), fardiag_val));
		coefficients.push_back(Trip(idx + (N - 2), idx, gamma));
		coef2.push_back(Trip(idx + (N - 2), idx, fardiag_val));
	}


	// i = 1 // i = N - 2 // j = 1 // j = N - 2
	for (int j = 1; j < N - 1; j++) {
		int idx = to_idx(1, j, N);
		B[idx] += -phi_1(y_step * j) / (x_step * x_step);

		idx = to_idx(N - 2, j, N);
		B[idx] += -phi_2(y_step * j) / (x_step * x_step);

		idx = to_idx(j, 1, N);
		B[idx] += -phi_3(x_step * j) / (y_step * y_step);

		idx = to_idx(j, N - 2, N);
		B[idx] += -phi_4(x_step * j) / (y_step * y_step);
	}


	SpMat A((N - 2) * (N - 2), (N - 2) * (N - 2));
	A.setFromTriplets(coefficients.begin(), coefficients.end());

	SpMat D((N - 2) * (N - 2), (N - 2) * (N - 2));
	D.setFromTriplets(coef2.begin(), coef2.end());

	{
		// calculating ||I-DA||
		std::vector<Trip> I_tr;
		Eigen::VectorXd u_cur((N - 2) * (N - 2));
		for (int idx = 0; idx < (N - 2) * (N - 2); idx++) {
			I_tr.push_back(Trip(idx, idx, 1));
			u_cur[idx] = 0.;
		}
		SpMat I((N - 2) * (N - 2), (N - 2) * (N - 2));
		I.setFromTriplets(I_tr.begin(), I_tr.end());
		double koef = +0.0001;
		SpMat res = I - D * A;
		//Eigen::VectorXcd eig_vals;
		{
			/*
			auto _t = Eigen::MatrixXd(res);
			Eigen::EigenSolver<Eigen::MatrixXd> _es;
			_es.compute(_t, false);
			auto eig_vals = _es.eigenvalues();
			double max_eig_val = 0;
			for (int i = 0; i < eig_vals.size(); i++) {
				max_eig_val = max(abs(eig_vals[i].real()), max_eig_val);
			}

			cout << "max eig_val = " << max_eig_val << endl;
			*/
		}

		//for (int i = 0; i < eig_vals.size(); i++) {
			//cout << eig_vals[i] << endl;
		//}

		for (int i = 0; i < 50; i++) {
			u_cur = res * u_cur + koef * B;
		}
		ofstream f_out(out_filename + "_iters");
		for (int i = 1; i < N - 1; i++) {
			for (int j = 1; j < N - 1; j++) {
				f_out << u_cur[to_idx(i, j, N)] << " ";
			}
			f_out << endl;
		}
		f_out.close();
	}

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
	Eigen::VectorXd u_eig = chol.solve(B);        // use the factorization to solve for the given right hand side

	///////////////////////////////////////////


	// save res
	ofstream f_out(out_filename + "_pars");
	f_out << N << endl << iter_num << endl << err << endl;
	f_out.close();
	f_out.open(out_filename + "_true_vals");
	for (int i = 1; i < N - 1; i++) {
		for (int j = 1; j < N - 1; j++) {
			f_out << u_true[to_idx(i, j, N)] << " ";
		}
		f_out << endl;
	}
	f_out.close();
	f_out.open(out_filename + "_eig_vals");
	for (int i = 1; i < N - 1; i++) {
		for (int j = 1; j < N - 1; j++) {
			f_out << u_eig[to_idx(i, j, N)] << " ";
		}
		f_out << endl;
	}
	f_out.close();
}