#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct> 
#include <unsupported/Eigen/MatrixFunctions>
#include "GetRefPath.h"
#include "osqp/osqp.h"
#include "OsqpEigen/OsqpEigen.h"
#include <memory>

using namespace Eigen;

class MPCCal {
public:
	MPCCal(VectorXd sta_local);
	void MPCInit();
	bool MPCSolver(GetRefPath* ref_ptr);
	void ConsCal(SparseMatrix<double>& Acon_sparse, VectorXd& ub);
	vector<VectorXd> status_node;
private:
	double safe_r = 0.25; //°²È«·¶Î§
	double obs_r = 0.5; //ÕÏ°­Îï°ë¾¶

	vector<Vector2d> obs_vec; 

	Matrix<double, STANUM, STANUM> A;
	Matrix<double, STANUM, CONNUM> B;

	VectorXd sta_local = VectorXd::Zero(STANUM,1);

	MatrixXd Q = MatrixXd::Zero(STANUM, STANUM);
	MatrixXd Q_bar = MatrixXd::Zero(STANUM * N, STANUM * N);
	MatrixXd R = MatrixXd::Identity(CONNUM, CONNUM);
	MatrixXd R_bar = MatrixXd::Zero(CONNUM * N, CONNUM * N);

	VectorXd u = VectorXd::Zero(CONNUM);
	VectorXd solve = VectorXd::Zero(CONNUM*N);

	MatrixXd M = MatrixXd::Zero(N * STANUM, STANUM);
	MatrixXd T = MatrixXd::Zero(N * STANUM, N * CONNUM);

	MatrixXd H = MatrixXd::Zero(N * CONNUM, N * CONNUM);
	SparseMatrix<double> H_sparse;

	VectorXd lb = VectorXd::Constant(CONSROWNUM, -OSQP_INFTY);

	Vector2d obs1 = {4.0,4.0};
};