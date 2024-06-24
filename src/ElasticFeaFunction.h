#pragma once
#ifndef ELASTICFEAFUNCTION_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define ELASTICFEAFUNCTION_H

#include "FeaHelperFunctions.h"

Eigen::MatrixXd ElasticFea(const Eigen::MatrixXd Rho, const MESH myMesh, const Eigen::MatrixXd Hmat, const Eigen::MatrixXd T, FEA& myFea, Eigen::SparseMatrix<double>& Kg, GRAD_HELPER_VARS& myGrad_Helper_Vars);

#endif