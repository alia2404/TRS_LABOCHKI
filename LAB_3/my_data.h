#pragma once

#include <cmath>
#include <vector>
#define PI acos(-1.)

#define LX 1.
#define LY 2.

#define ITER_CONST 0
#define ERR_CONST 1
using namespace std;

static double u_func(double x, double y) {
	return sin(PI * x * x) * sin(PI * y / 2);
}

static double phi_1(double y) {
	return 0; // u_func(0, y);
}

static double phi_2(double y) {
	return 0; // u_func(LX, y);
}

static double phi_3(double x) {
	return 0; // u_func(x, 0);
}

static double phi_4(double x) {
	return 0; // u_func(x, LY);
}

static double f_func_1(double x, double y) {
	return PI * PI * (1 + 16 * x * x) * sin(PI * x * x) * sin(PI * y / 2) / 4
		- 2 * PI * cos(PI * x * x) * sin(PI * y / 2);
}

static double a_func(double x, double y) {
	return 1;
}

static double b_func(double x, double y) {
	return 1 + y * y;
}

static double c_func(double x, double y) {
	return x * x - y * y;
}

static double f_func_2(double x, double y) {
	return
		u_func(x, y) * (4 * PI * PI * x * x + PI * PI / 4 + PI * PI * y * y / 4 - c_func(x, y))
		- 2 * PI * cos(PI * x * x) * sin(PI * y / 2)
		- PI * y * sin(PI * x * x) * cos(PI * y / 2);
}

static int to_idx(int i, int j, int N) {
	return (j - 1) * (N - 1) + (i - 1);
}

static void to_ij(int idx, int N, int& i, int& j) {
	j = idx / (N - 1) + 1;
	i = idx % (N - 1) + 1;
}

static double calc_delta(std::vector<double> v1, std::vector<double> v2) {
	double delta = abs(v1[0] - v2[0]);
	for (int i = 1; i < v2.size(); i++) {
		delta = std::max(delta, abs(v1[i] - v2[i]));
	}
	return delta;
}