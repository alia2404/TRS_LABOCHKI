#pragma once

#include "eigen-3.4.0\Eigen\Sparse"

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Trip;

void t3_eigen(double x_step, double y_step, double& err, int& iter_num, int mode, std::string out_filename);
void t3(double x_step, double y_step, double& err, int& iter_num, int mode, std::string out_filename);