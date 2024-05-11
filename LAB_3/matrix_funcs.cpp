#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <stdlib.h> 
#include "matrix_funcs.h"
#include "my_data.h"

#define ITER_LIM 100

using namespace std;

#pragma region Вспомогательные функции

vector<vector<double>> random_matrix(int M, int N) {
	vector<vector<double>> A(M, vector<double>(N, 0));
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			A[i][j] = 2 * ((double)rand()) / RAND_MAX - 1;
		}
	}
	return A;
}

// Произведение матрицы на число
vector<vector<double>> matrix_prod(vector <vector <double>> A, double C) {
	vector<vector<double>> A1 = A;
	for (int j = 0; j < A.size(); j++)
		for (int k = 0; k < A[0].size(); k++)
			A1[j][k] *= C;
	return A1;
}

// Произведение матриц
vector<vector<double>> matrix_prod(vector <vector <double>> A, vector <vector <double>> B)
{
	vector<vector<double>> R;
	if (A[0].size() != B.size()) return R;
	R.resize(A.size());
	for (int i = 0; i < A.size(); i++) {
		R[i].resize(B[0].size(), 0);
		for (int j = 0; j < B[0].size(); j++) {
			for (int k = 0; k < B.size(); k++) {
				R[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	return R;
}

// Сумма матриц
vector<vector<double>> matrix_add(vector <vector <double>> A, vector <vector <double>> B)
{
	vector<vector<double>> R;
	if (A.size() != B.size() || A[0].size() != B[0].size()) return R;
	R.resize(A.size());
	for (int i = 0; i < A.size(); i++) {
		R[i].resize(A[0].size(), 0);
		for (int j = 0; j < A[0].size(); j++) {
			R[i][j] = A[i][j] + B[i][j];
		}
	}
	return R;
}

// Разность матриц
vector<vector<double>> matrix_sub(vector <vector <double>> A, vector <vector <double>> B)
{
	vector<vector<double>> R;
	if (A.size() != B.size() || A[0].size() != B[0].size()) return R;
	R.resize(A.size());
	for (int i = 0; i < A.size(); i++) {
		R[i].resize(A[0].size(), 0);
		for (int j = 0; j < A[0].size(); j++) {
			R[i][j] = A[i][j] - B[i][j];
		}
	}
	return R;
}


// Обратная к нижнетреугольной
vector<vector<double>> low_tr_matrix_inv(vector <vector <double>> A) {
	int N = A.size();
	vector<vector<double>> A_inv(N, vector<double>(N, 0));
	for (int i = 0; i < N; i++) {
		A_inv[i][i] = 1;
	}
	for (int k = 0; k < N - 1; k++) {
		for (int i = k + 1; i < N; i++) {
			for (int j = 0; j < k + 1; j++) {
				A_inv[i][j] -= A_inv[k][j] * A[i][k] / A[k][k];
			}
		}
		for (int j = 0; j < k + 1; j++) {
			A_inv[k][j] /= A[k][k];
		}
	}
	for (int j = 0; j < N; j++) {
		A_inv[N - 1][j] /= A[N - 1][N - 1];
	}
	return A_inv;
}

// Обратная к верхнетреугольной
vector<vector<double>> up_tr_matrix_inv(vector <vector <double>> A) {
	int N = A.size();
	vector<vector<double>> A_inv(N, vector<double>(N, 0));
	for (int i = 0; i < N; i++) {
		A_inv[i][i] = 1;
	}
	for (int k = N - 1; k > 0; k--) {
		for (int i = k - 1; i > -1; i--) {
			for (int j = k; j < N; j++) {
				A_inv[i][j] -= A_inv[k][j] * A[i][k] / A[k][k];
			}
		}
		for (int j = k; j < N; j++) {
			A_inv[k][j] /= A[k][k];
		}
	}
	for (int j = 0; j < N; j++) {
		A_inv[0][j] /= A[0][0];
	}
	return A_inv;
}

// Обратное вычисление снизу вверх
void back_substitution(vector<vector<double>> A, vector<double>& X, vector<double> B) {
	X.resize(B.size(), 0);
	for (int i = X.size() - 1; i >= 0; i--) {
		double sum = 0.;
		for (int j = i + 1; j < X.size(); j++) {
			sum += A[i][j] * X[j];
		}
		X[i] = (B[i] - sum) / A[i][i];
	}
}


// Разделяет матрицу на три части
void separate_matrix(vector <vector <double>> A,
	vector <vector <double>>& D,
	vector <vector <double>>& E,
	vector <vector <double>>& F) {
	cout << "start separate" << endl;
	int N = A.size();
	cout << "1 separate" << endl;
	D.resize(N, vector<double>(N, 0));
	cout << "2 separate" << endl;
	E.resize(N, vector<double>(N, 0));
	cout << "3 separate" << endl;
	F.resize(N, vector<double>(N, 0));
	cout << "4 separate" << endl;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < i; j++) {
			E[i][j] = -A[i][j];
		}
		D[i][i] = A[i][i];
		for (int j = i + 1; j < N; j++) {
			F[i][j] = -A[i][j];
		}
	}
}


double scal_prod(vector<double> x, vector<double> y) {
	double prod = 0;
	for (int i = 0; i < x.size(); i++) {
		prod += x[i] * y[i];
	}
	return prod;
}

double scal_prod(vector< vector<double>> x,
	vector< vector<double>> y) {
	if (x[0].size() == 1) {
		double prod = 0;
		for (int i = 0; i < x.size(); i++) {
			prod += x[i][0] * y[i][0];
		}
		return prod;
	}
	else {
		double prod = 0;
		for (int i = 0; i < x[0].size(); i++) {
			prod += x[0][i] * y[0][i];
		}
		return prod;
	}
}

double eucl_norm(vector<double> x) {
	return sqrt(scal_prod(x, x));
}

double eucl_norm(vector<vector<double>> x) {
	double norm = 0;
	for (int i = 0; i < x.size(); i++) {
		norm += pow(eucl_norm(x[i]), 2.);
	}
	return sqrt(norm);
}


// Транспонирование
void matrix_transpose(vector <vector <double>> A, vector <vector <double>>& A_T) {
	A_T.clear();
	A_T.resize(A[0].size());
	for (int i = 0; i < A_T.size(); i++) {
		A_T[i].resize(A.size());
	}

	for (int i = 0; i < A.size(); i++) {
		for (int j = 0; j < A[0].size(); j++) {
			A_T[j][i] = A[i][j];
		}
	}
}

// Обратное вычисление сверху вниз
void forward_substitution(vector<vector<double>> A, vector<double>& X, vector<double> B) {
	for (int i = 0; i < A.size(); i++) {
		double sum = 0.;
		for (int j = 0; j < i; j++) {
			sum += A[i][j] * X[j];
		}
		X[i] = (B[i] - sum) / A[i][i];
	}
}

#pragma endregion

#pragma region SOR

vector<vector<double>> sor_step(vector <vector <double>> X,
	vector <vector <double>> D,
	vector <vector <double>> F,
	vector <vector <double>> B,
	vector <vector <double>> DmwE_inv,
	double w) {
	auto G_w = matrix_prod(
		DmwE_inv,
		matrix_add(matrix_prod(F, w), matrix_prod(D, 1 - w)));

	auto f_w = matrix_prod(
		matrix_prod(DmwE_inv, B),
		w);

	return matrix_add(
		matrix_prod(G_w, X),
		f_w
	);
}

vector <vector <double>> solve_sor(vector <vector <double>> A,
	vector <vector <double>> B,
	double err, double w,
	int& iter_number) {
	int N = A.size();
	vector <vector <double>> X = random_matrix(N, 1);
	vector <vector <double>> D;
	vector <vector <double>> E;
	vector <vector <double>> F;
	separate_matrix(A, D, E, F);

	auto DmwE_inv = low_tr_matrix_inv(matrix_sub(D, matrix_prod(E, w)));
	auto X_new = sor_step(X, D, F, B, DmwE_inv, w);
	auto delta_X = matrix_sub(X, X_new);
	double delta = eucl_norm(delta_X);
	iter_number = 1;
	while (delta > err && iter_number <= ITER_LIM) {
		X = X_new;
		X_new = sor_step(X, D, F, B, DmwE_inv, w);
		delta_X = matrix_sub(X, X_new);
		delta = eucl_norm(delta_X);
		iter_number++;
	}
	return X_new;
}



// SOR. Известно ист. значение.
vector <vector <double>> solve_sor(vector <vector <double>> A,
	vector <vector <double>> B,
	double err, double w,
	vector <vector <double>> X_true, int& iter_number) {
	int N = A.size();
	vector <vector <double>> X = random_matrix(N, 1);
	vector <vector <double>> D;
	vector <vector <double>> E;
	vector <vector <double>> F;
	separate_matrix(A, D, E, F);
	auto DmwE_inv = low_tr_matrix_inv(matrix_sub(D, matrix_prod(E, w)));
	auto X_new = sor_step(X, D, F, B, DmwE_inv, w);
	auto delta_X = matrix_sub(X, X_new);
	double delta = eucl_norm(delta_X);
	iter_number = 1;
	while (delta > err && iter_number <= ITER_LIM) {
		X = X_new;
		X_new = sor_step(X, D, F, B, DmwE_inv, w);
		delta_X = matrix_sub(X_new, X_true);
		delta = eucl_norm(delta_X);
		iter_number++;
	}
	return X_new;
}

vector<vector<double>> sym_sor_step(vector <vector <double>> X,
	vector <vector <double>> D,
	vector <vector <double>> E,
	vector <vector <double>> F,
	vector <vector <double>> B,
	vector <vector <double>> DmwE_inv,
	vector <vector <double>> DmwF_inv,
	double w) {
	auto G_w = matrix_prod(
		matrix_prod(
			DmwF_inv,
			matrix_add(matrix_prod(E, w), matrix_prod(D, 1 - w))),
		matrix_prod(
			DmwE_inv,
			matrix_add(matrix_prod(F, w), matrix_prod(D, 1 - w)))
	);

	auto f_w = matrix_prod(
		matrix_prod(
			matrix_prod(DmwF_inv, D),
			matrix_prod(DmwE_inv, B)),
		w * (2 - w));

	return matrix_add(
		matrix_prod(G_w, X),
		f_w
	);
}

// SOR. Неизвестно ист. значение.
vector <vector <double>> solve_sym_sor(vector <vector <double>> A,
	vector <vector <double>> B,
	double& err, double w,
	int& iter_num, int mode) {
	int N = A.size();
	vector <vector <double>> X = random_matrix(N, 1);
	vector <vector <double>> D(N);
	vector <vector <double>> E(N);
	vector <vector <double>> F(N);
	//separate_matrix(A, D, E, F);
	
	{
		for (int i = 0; i < N; i++) {
			//cout << "i=" << i << endl;
			E[i] = vector<double>(N, 0);
			D[i] = vector<double>(N, 0);
			F[i] = vector<double>(N, 0);
			for (int j = 0; j < i; j++) {
				E[i][j] = -A[i][j];
			}
			D[i][i] = A[i][i];
			for (int j = i + 1; j < N; j++) {
				F[i][j] = -A[i][j];
			}
		}
	}

	auto DmwE_inv = low_tr_matrix_inv(matrix_sub(D, matrix_prod(E, w)));
	auto DmwF_inv = up_tr_matrix_inv(matrix_sub(D, matrix_prod(F, w)));
	auto X_new = X;
	int iter_cnt = 0;
	double curr_err = 0;
	do {
		if ((iter_cnt+1) % 1 == 0) cout << "Iter number " << iter_cnt + 1 << " START" << endl;
		X = X_new;
		X_new = sym_sor_step(X, D, E, F, B, DmwE_inv, DmwF_inv, w);

		iter_cnt++;
		if (mode == ITER_CONST) {
			if (iter_cnt >= iter_num) {
				auto delta_X = matrix_sub(X, X_new);
				err = eucl_norm(delta_X);
				break;
			}
		}
		else /*if (mode == ERR_CONST)*/ {
			auto delta_X = matrix_sub(X, X_new);
			curr_err = eucl_norm(delta_X);
			if (curr_err <= err) {
				iter_num = iter_cnt;
				break;
			}
		}
		if (iter_cnt % 1 == 0) cout << "Iter number " << iter_cnt << "; max_diff = " << curr_err << endl;
	} while (true);
	return X_new;
}

// SOR. Известно ист. значение.
vector <vector <double>> solve_sym_sor(vector <vector <double>> A,
	vector <vector <double>> B,
	double err, double w,
	vector <vector <double>> X_true, int& iter_number) {
	int N = A.size();
	vector <vector <double>> X = random_matrix(N, 1);
	vector <vector <double>> D;
	vector <vector <double>> E;
	vector <vector <double>> F;
	separate_matrix(A, D, E, F);
	auto DmwE_inv = low_tr_matrix_inv(matrix_sub(D, matrix_prod(E, w)));
	auto DmwF_inv = up_tr_matrix_inv(matrix_sub(D, matrix_prod(F, w)));
	auto X_new = sym_sor_step(X, D, E, F, B, DmwE_inv, DmwF_inv, w);
	auto delta_X = matrix_sub(X_new, X_true);
	double delta = eucl_norm(delta_X);
	iter_number = 1;
	while (delta > err) {
		X = X_new;
		X_new = sym_sor_step(X, D, E, F, B, DmwE_inv, DmwF_inv, w);
		delta_X = matrix_sub(X_new, X_true);
		delta = eucl_norm(delta_X);
		iter_number++;
	}
	return X_new;
}
#pragma endregion

#pragma region Square_root_method
void s_matrix_define(vector<vector<double>> A, vector<vector<double>>& S, vector<vector<double>>& SIGN) {
	double sum = 0., bet = 0.;

	SIGN[0][0] = A[0][0] > 0 ? 1 : -1;
	S[0][0] = pow(abs(A[0][0]), 0.5);

	for (int k = 0; k < A.size(); k++) {
		sum = 0;

		for (int i = 0; i < k; i++) {
			sum += SIGN[i][i] * S[i][k] * S[i][k];
		}
		bet = A[k][k] - sum;

		SIGN[k][k] = bet > 0 ? 1 : -1;
		S[k][k] = pow(abs(bet), 0.5);

		for (int l = k + 1; l < A.size(); l++) {
			sum = 0;
			for (int i = 0; i < k; i++) {
				sum += SIGN[i][i] * S[i][k] * S[i][l];
			}

			S[k][l] = (A[k][l] - sum) * SIGN[k][k] / (S[k][k]);
		}
	}
}

void Square_root_method(vector<vector<double>> A, vector<double> B, vector<double>& X) {
	vector<vector<double>> S(A.size(), vector<double>(A.size())),
		S_T(A.size(), vector<double>(A.size()));
	vector<vector<double>> SIGN(A.size(), vector<double>(A.size()));
	vector<double> Y(A.size());

	s_matrix_define(A, S, SIGN);
	matrix_transpose(S, S_T);

	forward_substitution(S_T, Y, B);
	back_substitution(S, X, Y);

	//show_X(X);
}
#pragma endregion

#pragma region LU_decay
void lu_matrix(vector <vector <double>> A, vector <vector <double>>& L, vector <vector <double>>& U) {
	U = A;
	for (int i = 0; i < A.size(); i++)
		for (int j = i; j < A.size(); j++)
			L[j][i] = U[j][i] / U[i][i];


	for (int k = 1; k < A.size(); k++) {
		for (int i = k - 1; i < A.size(); i++)
			for (int j = i; j < A.size(); j++)
				L[j][i] = U[j][i] / U[i][i];

		for (int i = k; i < A.size(); i++)
			for (int j = k - 1; j < A.size(); j++)
				U[i][j] = U[i][j] - L[i][k - 1] * U[k - 1][j];
	}
}

void LU_decay(vector<vector<double>> A, vector<double> B, vector<double>& X) {
	vector<vector<double>> L(A.size(), vector<double>(A.size())),
		U(A.size(), vector<double>(A.size()));
	vector<double> Y(A.size());

	lu_matrix(A, L, U);

	forward_substitution(L, Y, B);
	back_substitution(U, X, Y);

	//show_X(X);
}
#pragma endregion