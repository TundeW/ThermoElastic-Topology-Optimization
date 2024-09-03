#pragma once
#ifndef FOWADANALYSISFUNCTIONS_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define FOWARDANALYSISFUNCTIONS_H

//#include <Eigen/Dense>
#include "ProjectionHelperFunctions.h"

typedef struct {
    double MaxStress;
	double HeatTransferRate;
	double PressureDrop;
	double HF_VolFrac;
	Eigen::MatrixXd dMaxStress_dx;
	Eigen::MatrixXd dHeatTransferRate_dx;
	Eigen::MatrixXd dPressureDrop_dx;
	Eigen::MatrixXd circDom;
	Eigen::MatrixXd dcircDom_dx;
	Eigen::MatrixXd dHF_VolFracdx;

} CAND_OBJ_GRAD;

double ComputeAll(const Eigen::VectorXd x, const MESH myMesh, CAND_OBJ_GRAD& myCAND_OBJ_GRAD); //, int el, double h);


#endif