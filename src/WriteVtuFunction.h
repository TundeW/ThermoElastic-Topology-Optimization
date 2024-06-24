#pragma once
#ifndef WRITEVTUFUNCTION_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define WRITEVTUFUNCTION_H

#include "FeaHelperFunctions.h"
#include "ProjectionHelperFunctions.h"

void WriteVtu(const Eigen::MatrixXd Rho, const MESH myMesh, const FEA myFea);

#endif