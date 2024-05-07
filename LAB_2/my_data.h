#pragma once

#include <cmath>
#include <vector>
#define PI acos(-1.)

static double func(double x, double t) { return -x * x + x + 2. * t; }
static double phi(double x) { return sin(PI * x) / PI; }
static double psi0(double t) { return -exp(-PI * PI * t) - t; }
static double psi1(double t) { return -exp(-PI * PI * t) - t; }
static double Func(double u) { return 0+1*sin(u); }
static double k(double u) { return 0+1*u * u; }
static double sol(double x, double t) {
    return exp(-PI * PI * t) * sin(PI * x) / PI + t * x * (1 - x);
}

static void diag3_method(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d, int n, std::vector<double>& x) {
    std::vector<double> ksi = std::vector<double>(n + 1);
    std::vector<double> teta = std::vector<double>(n + 1);
    ksi[0] = 0;
    teta[0] = 0;
    for (int i = 0; i < n; i++) {
        double denominator = b[i] - a[i] * ksi[i];
        ksi[i + 1] = c[i] / denominator;
        teta[i + 1] = (a[i] * teta[i] - d[i]) / denominator;
    }
    x = std::vector<double>(n);
    x[n - 1] = teta[n];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = ksi[i + 1] * x[i + 1] + teta[i + 1];
    }
}
