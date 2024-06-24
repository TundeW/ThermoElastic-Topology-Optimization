#pragma once
#ifndef THERMALFEAFUNCTION_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define THERMALFEAFUNCTION_H

#include "FeaHelperFunctions.h"

// Eigen::MatrixXd ThermalFea(const Eigen::MatrixXd Rho, const MESH myMesh, const Eigen::MatrixXd Hmat);
Eigen::MatrixXd ThermalFea(const Eigen::MatrixXd Rho, const MESH myMesh, const Eigen::SparseMatrix<double> Hmat, FEA& myFea, Eigen::SparseMatrix<double>& Kg, GRAD_HELPER_VARS& myGrad_Helper_Vars);

#endif