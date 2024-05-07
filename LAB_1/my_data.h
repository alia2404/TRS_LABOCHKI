#pragma once
#include <cmath>

static double m = 1.67e-27;
static double x0 = 0;
static double v0 = 7e5;
static double a = 6.396483410;
static double b = 700000.;
static double c = 9.137833443e-6;

static double u(double x) { return 2e-17 * cosh(x); }
static double du(double x) { return 2e-17 * sinh(x); }
static double x_sc(double x) { return x / a; } // sc = scale
static double v_sc(double v) { return v / b; }
static double t_sc(double t) { return t / c; }
static double x_or(double x_sc) { return x_sc * a; } // or = origin
static double v_or(double v_sc) { return v_sc * b; }
static double t_or(double t_sc) { return t_sc * c; }
static double dxdt_sc(double x_sc, double v_sc) { return b * c * v_sc / a; }
static double dvdt_sc(double x_sc, double v_sc) { return -c * 2e-17 * sinh(a * x_sc) / (b * m); }
