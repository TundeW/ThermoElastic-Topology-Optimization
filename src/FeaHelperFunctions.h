#pragma once
#ifndef FEAHELPERFUNCTIONS_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define FEAHELPERFUNCTIONS_H

using namespace Eigen;
#include "ProjectionHelperFunctions.h"

typedef struct {
    Eigen::MatrixXd T;
	Eigen::MatrixXd U;
	Eigen::MatrixXd D_Elastic;
	Eigen::MatrixXd strain_thermal;
	Eigen::MatrixXd vonMises;
	Eigen::MatrixXd stress_true;
    Eigen::MatrixXd Vel;
	Eigen::MatrixXd Kappa_Tau;
	double K1;
	/*double nel;
	int dim;
	double L;
	double H;
	double B;
	double nelx;
	double nely;
	double nelz;
	double Area;*/
} FEA;

void shapeFcn(const int i, const Eigen::MatrixXd gPoints, const MESH myMesh, Eigen::MatrixXd& N, Eigen::MatrixXd& dN);
//void elem_stiff(const Eigen::MatrixXd D, const Eigen::MatrixXd elem_coords, const Eigen::MatrixXd gPoints, const Eigen::MatrixXd strain_thermal, const Eigen::MatrixXd wt, const double A, const MESH myMesh, Eigen::MatrixXd& k_el, Eigen::MatrixXd& f_thermal_el);
void elem_stiff(const Eigen::MatrixXd D, const Eigen::MatrixXd elem_coords, const Eigen::MatrixXd gPoints, const Eigen::MatrixXd strain_thermal, const Eigen::MatrixXd wt, const double A, const MESH myMesh, MatrixXd& k_el, MatrixXd& f_thermal_el, MatrixXd& df_thermal_eldt, MatrixXd& df_thermal_eldstrain_thermal, const Eigen::MatrixXd dstrain_thermaldt, MatrixXd& dk_el_drho1, MatrixXd& dk_el_drho2, const Eigen::MatrixXd dD_rho_drho1, const Eigen::MatrixXd dD_rho_drho2, MatrixXd& delf_thermal_el_delrho1, MatrixXd& delf_thermal_el_delrho2);
Eigen::MatrixXd nodes_to_dofs(const Eigen::MatrixXd enodes, const MESH myMesh);
double ComputeStress(const Eigen::MatrixXd Rho, const MESH myMesh, const Eigen::MatrixXd Hmat, FEA& myFea, GRAD_HELPER_VARS& myGrad_Helper_Vars);
double ComputeHeatTransferRate(const Eigen::MatrixXd Rho, const MESH myMesh, FEA& myFea, GRAD_HELPER_VARS& myGrad_Helper_Vars);
//std::vector<std::pair<std::string, std::vector<int>>> read_csv(std::string filename);
std::vector<std::vector<double>> read_csv(const std::string& filename);

#endif