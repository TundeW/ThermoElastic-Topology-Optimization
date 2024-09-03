#pragma once
#ifndef PROJECTIONHELPERFUNCTIONS_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define PROJECTIONHELPERFUNCTIONS_H

#include <Eigen/Dense>

struct MESH {
    Eigen::MatrixXd coordinate;
    Eigen::MatrixXd nodemap;
    Eigen::MatrixXd mesh_centers;
    double nel;
    int dim;
    int conLen;
    int opt_dim;
    int opt_iter;
    double L;
    double H;
    double B;
    double nelx;
    double nely;
    double nelz;
    double Area;
    double he;
    Eigen::SparseMatrix<double> Hmat;
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
    
};

typedef struct {
    Eigen::MatrixXd dRho1_dx;
	Eigen::MatrixXd dRho2_dx;
	std::vector<Eigen::MatrixXd> dKthdrho1;
	std::vector<Eigen::MatrixXd> dKthdrho2;
	std::vector<Eigen::MatrixXd> dPthdrho1;
	std::vector<Eigen::MatrixXd> dPthdrho2;
    std::vector<Eigen::MatrixXd> dKeldrho1;
	std::vector<Eigen::MatrixXd> dKeldrho2;
	std::vector<Eigen::MatrixXd> dPeldrho1;
	std::vector<Eigen::MatrixXd> dPeldrho2;
	std::vector<Eigen::MatrixXd> dstrain_thermaldrho1;
	std::vector<Eigen::MatrixXd> dstrain_thermaldrho2;
	std::vector<Eigen::MatrixXd> dstrain_thermaldt;
    Eigen::SparseMatrix<double> dFdT;
    Eigen::SparseMatrix<double> Kth_freefree;
    Eigen::MatrixXi th_freedofs;
    Eigen::SparseMatrix<double> Kel_freefree;
    Eigen::MatrixXi el_freedofs;
    Eigen::MatrixXd dfdrho1;
    Eigen::MatrixXd dfdrho2;
    Eigen::MatrixXd df2drho1;
    Eigen::MatrixXd df2drho2;

} GRAD_HELPER_VARS;

Eigen::MatrixXd BarProjection(const Eigen::MatrixXd Centre, const Eigen::VectorXd R, const MESH myMesh, Eigen::MatrixXd& dRhoodx, Eigen::MatrixXd& dRhoodr, Eigen::MatrixXd& dRhood_);
Eigen::MatrixXd PipeProjection(const Eigen::MatrixXd Centre, const Eigen::VectorXd Inner_r, const Eigen::VectorXd Outer_r, const MESH myMesh, GRAD_HELPER_VARS& myGrad_Helper_Vars);

#endif