#pragma once

#include <Eigen\Sparse>

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Trip;

void t4(double x_step, double y_step, double& err, int& iter_num, int mode, std::string out_filename);