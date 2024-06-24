#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <iostream>
#include <fstream>
#include "FowardAnalysisFunctions.h"
#include "ProjectionHelperFunctions.h"
#include "FeaHelperFunctions.h"
#include "ThermalFeaFunction.h"
#include "ElasticFeaFunction.h"
#include "WriteVtuFunction.h"
#include <omp.h>
using namespace Eigen;

//double ComputeAll(const Eigen::VectorXd x, const MESH myMesh, CAND_OBJ_GRAD& myCAND_OBJ_GRAD, int el, double h) {
double ComputeAll(const Eigen::VectorXd x, const MESH myMesh, CAND_OBJ_GRAD& myCAND_OBJ_GRAD) {
    // Access struct parameters
	int dim = myMesh.dim;
	Eigen::SparseMatrix<double> Hmat = myMesh.Hmat;

	int pipenum = x.size()/4;
	// Define N and dN outside the if statement, using dim
	std::cout << "Defining Design Variables (Pipe Parameters): Location, Inner Radius and Pipe Thickness" << std::endl;
	/*MatrixXd Centre(pipenum, 2);
   	VectorXd Inner_r(pipenum);
	VectorXd Outer_r(pipenum);
   	Centre << 0.025, 0.025, 0.175, 0.025; //, 0.175, 0.025, 0.075, 0.075, 0.125, 0.075, 0.225, 0.075;
   	Inner_r << .01, .01; //, .01, .01, .01, .01; //, .01, .01, .01, .01;
	Outer_r << .02, .02; //, .02, .02, .02, .02; //, .02, .02, .02, .02;*/
	//0.125, 0.025, 0.225, 0.025, 0.025, 0.075, 0.175, 0.075,
	//Centre << 0.05, 0.05;
	//Inner_r << .02;
	//Outer_r << .04;

	
	MatrixXd centre = x.head(pipenum*2).reshaped(2, pipenum);
	MatrixXd Centre = centre.transpose();
    VectorXd Inner_r = x.segment(pipenum*2, pipenum);
    VectorXd Outer_r = x.tail(pipenum);

	//std::cout << "Centre:\n" << Centre << std::endl;
	//std::cout << "Inner_r: " << Inner_r << std::endl;
	//std::cout << "Outer_r: " << Outer_r << std::endl;


   	std::cout << "Beginning Geometric Projection" << std::endl;
	auto start_time = std::chrono::high_resolution_clock::now();
	GRAD_HELPER_VARS myGrad_Helper_Vars;
	Eigen::MatrixXd Rho = PipeProjection(Centre, Inner_r, Outer_r, myMesh, myGrad_Helper_Vars);
	//std::cout << "el: " << el << std::endl;

	//Eigen::MatrixXd Rho = Eigen::MatrixXd::Constant(Rhop.rows(), Rhop.cols(), 0.5);
	//Rho(el,1) += h;
	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "PipeProjection time: " << duration.count() << " seconds" << std::endl;
	//std::cout << "Rho:\n" << Rho << std::endl;
    //std::cout << "Rho:\n" << Rho.rows() << " x " << Rho.cols() << "\n";
	//std::cout << "dRho1_dx:\n" << myGrad_Helper_Vars.dRho1_dx << "\n";
	//std::cout << "dRho2_dx:\n" << myGrad_Helper_Vars.dRho2_dx << "\n";
	FEA myFea;

	//start_time = std::chrono::high_resolution_clock::now();
	Eigen::MatrixXd nodemap = myMesh.nodemap;
	double nel = myMesh.nel;
	int ndof = nodemap.maxCoeff() + 1;
	Eigen::SparseMatrix<double> Kth(ndof, ndof);
	Kth.reserve(VectorXi::Constant(ndof,27));

	myGrad_Helper_Vars.dKthdrho1.resize(nel, Eigen::MatrixXd::Zero(pow(2,dim), pow(2,dim)));
	myGrad_Helper_Vars.dKthdrho2.resize(nel, Eigen::MatrixXd::Zero(pow(2,dim), pow(2,dim)));
	myGrad_Helper_Vars.dPthdrho1.resize(nel, Eigen::MatrixXd::Zero(pow(2,dim), 1));
	myGrad_Helper_Vars.dPthdrho2.resize(nel, Eigen::MatrixXd::Zero(pow(2,dim), 1));

	Eigen::MatrixXd T = ThermalFea(Rho, myMesh, Hmat, myFea, Kth, myGrad_Helper_Vars);
	//end_time = std::chrono::high_resolution_clock::now();
	//duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    //std::cout << "ThermalFea time: " << duration.count() << " seconds" << std::endl;
    //std::cout << "T:\n" << T.transpose() << "\n\n";
	myFea.T = T;

	double f2 = ComputeHeatTransferRate(Rho, myMesh, myFea, myGrad_Helper_Vars);

	std::cout << "f2: " << f2 << std::endl;

	/*std::cout << "T: " << T.transpose() << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho1[0]:\n" << myGrad_Helper_Vars.dKthdrho1[0] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho1[2]:\n" << myGrad_Helper_Vars.dKthdrho1[2] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho1[14]:\n" << myGrad_Helper_Vars.dKthdrho1[14] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho2[0]:\n" << myGrad_Helper_Vars.dKthdrho2[0] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho2[2]:\n" << myGrad_Helper_Vars.dKthdrho2[2] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho2[14]:\n" << myGrad_Helper_Vars.dKthdrho2[14] << std::endl;*/

	//std::cout << "Got Here 1" << std::endl;
	// start_time = std::chrono::high_resolution_clock::now();
	ndof = dim * (nodemap.maxCoeff() + 1);
	//std::cout << "Got Here 2" << std::endl;
	Eigen::SparseMatrix<double> Kel(ndof, ndof);
	Kel.reserve(VectorXi::Constant(ndof,81));
	myGrad_Helper_Vars.dKeldrho1.resize(nel, Eigen::MatrixXd::Zero(pow(2,dim), pow(2,dim)));
	myGrad_Helper_Vars.dKeldrho2.resize(nel, Eigen::MatrixXd::Zero(pow(2,dim), pow(2,dim)));
	myGrad_Helper_Vars.dPeldrho1.resize(nel, Eigen::MatrixXd::Zero(pow(2,dim), 1));
	myGrad_Helper_Vars.dPeldrho2.resize(nel, Eigen::MatrixXd::Zero(pow(2,dim), 1));
	myGrad_Helper_Vars.dstrain_thermaldrho1.resize(nel, Eigen::MatrixXd::Zero(6, 1));
	myGrad_Helper_Vars.dstrain_thermaldrho2.resize(nel, Eigen::MatrixXd::Zero(6, 1));
	myGrad_Helper_Vars.dstrain_thermaldt.resize(nel, Eigen::MatrixXd::Zero(6, pow(2,dim)));
	myGrad_Helper_Vars.dFdT.resize(ndof, ndof/dim);
	//std::cout << "ndof: " << ndof << ", dim: " << dim << ", ndof/dim: " << ndof/dim << std::endl;
	//std::cout << "Got Here 3" << std::endl;
	myGrad_Helper_Vars.dFdT.reserve(VectorXi::Constant(ndof/dim,81)); //,ndof/dim ******* Switched rows and cols
	//std::cout << "Got Here 4" << std::endl;

	Eigen::MatrixXd U = ElasticFea(Rho, myMesh, Hmat, T, myFea, Kel, myGrad_Helper_Vars);
	// end_time = std::chrono::high_resolution_clock::now();
	// duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    // std::cout << "ElasticFea time: " << duration.count() << " seconds" << std::endl;
    //std::cout << "U:\n" << U.transpose() << "\n\n";
	myFea.U = U;

	/*std::cout << "U:\n" << U.transpose() << std::endl;
	std::cout << "myGrad_Helper_Vars.dKeldrho1[0]:\n" << myGrad_Helper_Vars.dKeldrho1[0] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKeldrho1[2]:\n" << myGrad_Helper_Vars.dKeldrho1[2] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKeldrho1[14]:\n" << myGrad_Helper_Vars.dKeldrho1[14] << std::endl;
	std::cout << "myGrad_Helper_Vars.dFdT(8,3):\n" << myGrad_Helper_Vars.dFdT.coeff(8,3) << std::endl;
	std::cout << "myGrad_Helper_Vars.dFdT(11,6):\n" << myGrad_Helper_Vars.dFdT.coeff(11,6) << std::endl;
	std::cout << "myGrad_Helper_Vars.dFdT(32,10):\n" << myGrad_Helper_Vars.dFdT.coeff(32,10) << std::endl;
	std::cout << "myGrad_Helper_Vars.dFdT(70,27):\n" << myGrad_Helper_Vars.dFdT.coeff(70,27) << std::endl;*/

	start_time = std::chrono::high_resolution_clock::now();
	double vonMises = ComputeStress(Rho, myMesh, Hmat, myFea, myGrad_Helper_Vars);
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "ComputeStress time: " << duration.count() << " seconds" << std::endl;
    //std::cout << "myGrad_Helper_Vars.dfdrho2: " << myGrad_Helper_Vars.dfdrho2 << std::endl;
	//std::cout << "myGrad_Helper_Vars.dfdrho1: " << myGrad_Helper_Vars.dfdrho1 << std::endl;
	//std::cout << "myGrad_Helper_Vars.df2drho2: " << myGrad_Helper_Vars.df2drho2 << std::endl;
	//std::cout << "myGrad_Helper_Vars.df2drho1: " << myGrad_Helper_Vars.df2drho1 << std::endl;
	//myFea.vonMises = vonMises;

	MatrixXd dfdx = myGrad_Helper_Vars.dfdrho1 * myGrad_Helper_Vars.dRho1_dx + myGrad_Helper_Vars.dfdrho2 * myGrad_Helper_Vars.dRho2_dx;
	MatrixXd df2dx = myGrad_Helper_Vars.df2drho1 * myGrad_Helper_Vars.dRho1_dx + myGrad_Helper_Vars.df2drho2 * myGrad_Helper_Vars.dRho2_dx;

	/*myMesh.MaxStress = vonMises;
	myMesh.HeatTransferRate = f2;
	myMesh.dMaxStress_dx = dfdx;
	myMesh.dHeatTransferRate_dx = df2dx;*/

	myCAND_OBJ_GRAD.MaxStress = vonMises;
	myCAND_OBJ_GRAD.HeatTransferRate = f2;
	myCAND_OBJ_GRAD.dMaxStress_dx = dfdx;
	myCAND_OBJ_GRAD.dHeatTransferRate_dx = df2dx;

	std::cout << "dfdx: " << dfdx<< std::endl;

	WriteVtu(Rho,  myMesh, myFea);
	//Eigen::shutdownParallel(); 

	return vonMises;
	//return myFea.K1;
}

