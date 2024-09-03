#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <nlopt.hpp>
//#include <cmath>
#include "MeshHelperFunctions.h"
#include "FowardAnalysisFunctions.h"
#include "ProjectionHelperFunctions.h"
#include "FeaHelperFunctions.h"
#include "ThermalFeaFunction.h"
#include "ElasticFeaFunction.h"
#include "WriteVtuFunction.h"
#include <omp.h>

using namespace Eigen;

double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
	MESH *myMesh = reinterpret_cast<MESH*>(data);
	CAND_OBJ_GRAD myCAND_OBJ_GRAD;
	Eigen::Map<const Eigen::VectorXd> x_eigen(x.data(), x.size());

    double Rho = ComputeAll(x_eigen, *myMesh, myCAND_OBJ_GRAD); 

	myMesh->MaxStress = myCAND_OBJ_GRAD.MaxStress;
	myMesh->HeatTransferRate = myCAND_OBJ_GRAD.HeatTransferRate;
	myMesh->dMaxStress_dx = myCAND_OBJ_GRAD.dMaxStress_dx;
	myMesh->dHeatTransferRate_dx = myCAND_OBJ_GRAD.dHeatTransferRate_dx;

	Eigen::MatrixXd dMaxStress_dx = myMesh->dMaxStress_dx;
	std::cout << "Objectives: " << myMesh->MaxStress << std::endl;
	//std::cout << "Objectives gradient: \n" << std::endl;
    if (!grad.empty()) {
		for (int i = 0; i < x.size(); i++) {
			
			grad[i] = dMaxStress_dx(i);
			//std::cout << grad[i] << std::endl;
		} 	
    }
	
	
    return myMesh->MaxStress;

}

double myvconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
  	MESH *myMesh = reinterpret_cast<MESH*>(data);
	std::cout << "Constraint: " << myMesh->HeatTransferRate << std::endl;
	//std::cout << "Constraint gradient: \n" << std::endl;
  	if (!grad.empty()) {
		for (int i = 0; i < x.size(); i++) {
			grad[i] = myMesh->dHeatTransferRate_dx(i);
			//std::cout << grad[i] << std::endl;
		}
  	}
	
	
  	return myMesh->HeatTransferRate;
}

int main() {
	std::cout << "Thermo Elastic Solver" << std::endl;

	//omp_set_num_threads(4);
	#pragma omp parallel
	{
		std::cout << "This is thread: " << omp_get_thread_num() << "\n";
	}
	Eigen::initParallel(); // Initialize Eigen's parallel module

	// Creating Mesh
	double L = 0.25; //x
	//L = 1;
	double H = 0.1; //y
	//H = 4.5;
	double B = 0.256; //z
	//B = 0.5;
	double nelx = 50; double nely = 50; double nelz = 20;
	nelx = 8; nely = 6; nelz = 4;
	//nelx =40; nely = 20; nelz = 40;
	//nelx = 50; nely = 25; nelz = 25;
	//nelx =40; nely = 20; nelz = 20;
	//nelx = 25; nely = 15; nelz = 10;

	int dim = 3;

	std::cout << "Generating a " << dim << "D rectangular mesh" << std::endl;

	double dx = L/nelx;
	double dy = H/nely;
	double dz = B/nelz;
	double he = pow(pow(dx,2) + pow(dy,2) + pow(dz,2),0.5);
	//he = dx;
	double Area = dx * dy * dz;
	double nel = nelx*nely*nelz;

	auto start_time = std::chrono::high_resolution_clock::now();
	MatrixXd coordinates = generateGridCoordinates(L, H, B, nelx, nely, nelz);
	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "generateGridCoordinate time: " << duration.count() << " seconds" << std::endl;

	//std::cout << "Coordinates: " << coordinates << std::endl;
	//std::cout << "Coordinates:\n";
	//std::cout << coordinates.rows() << ", " << coordinates.cols() << std::endl;
	start_time = std::chrono::high_resolution_clock::now();
	Eigen::MatrixXd nodemap = generateNodeMap(nelx, nely, nelz);
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "generateNodeMap time: " << duration.count() << " seconds" << std::endl;
    //std::cout << "Nodemap:\n" << nodemap.rows() << ", " << nodemap.cols() << "\n";
	//std::cout << "Nodemap: " << nodemap << std::endl;

	start_time = std::chrono::high_resolution_clock::now();
	Eigen::MatrixXd mesh_centers = generateMeshCenters(nodemap, coordinates, nel);
    end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "generateMeshCenters time: " << duration.count() << " seconds" << std::endl;
    //std::cout << "Mesh Centers:\n" << mesh_centers.rows() << ", " << mesh_centers.cols() << "\n";
	//std::cout << "mesh_centers: " << mesh_centers << std::endl;

	start_time = std::chrono::high_resolution_clock::now();
	//Eigen::MatrixXd Hmat = generateElemToNodalDensityMatrix(nodemap, nel);
	Eigen::SparseMatrix<double> Hmat = generateElemToNodalDensityMatrix(nodemap, nel);
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "generateElemToNodalDensityMatrix time: " << duration.count() << " seconds" << std::endl;
    // std::cout << "Hmat:\n" << Hmat << "\n";
	//std::cout << "Hmat:\n" << Hmat.rows() << " x " << Hmat.cols() << "\n";

	MESH myMesh = {coordinates, nodemap, mesh_centers, nel, dim, L, H, B, nelx, nely, nelz, Area, he, Hmat};
	CAND_OBJ_GRAD myCAND_OBJ_GRAD;

	Eigen::VectorXd x(12);
    x << 0.025, 0.025, 0.025, 0.175, 0.175, 0.025, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02;
	//x << 0.025, 0.025, 0.01, 0.02;
	double h = 0;
	int el = 0;
	double Rho = ComputeAll(x, myMesh, myCAND_OBJ_GRAD); //, el, h);

	/*int opt_dim = 12; //Uncomment Here for FD
	std::vector<double> x_init(opt_dim);
  	x_init[0] = 0.025; x_init[1] = 0.025; x_init[2] = 0.025; x_init[3] = 0.175; x_init[4] = 0.175; x_init[5] = 0.025;
	x_init[6] = 0.01; x_init[7] = 0.01; x_init[8] = 0.01; x_init[9] = 0.02; x_init[10] = 0.02; x_init[11] = 0.02;

	
	std::vector<double> grad(opt_dim);
  	grad[0] = 0.025; grad[1] = 0.025; grad[2] = 0.025; grad[3] = 0.175; grad[4] = 0.175; grad[5] = 0.025;
	grad[6] = 0.01; grad[7] = 0.01; grad[8] = 0.01; grad[9] = 0.02; grad[10] = 0.02; grad[11] = 0.02;

	//double output = myvfunc(x_init, grad, myMesh);

	Eigen::MatrixXd df_dx;
	df_dx = MatrixXd::Zero(1, x.size());
	//df_dx = MatrixXd::Zero(1, nel);
	
	h = 0.77e-4;
	for (int i = 0; i < x.size(); i++) {
	//for (int i = 0; i < nel; i++) {
		Eigen::VectorXd x_pert = x;
		x_pert(i) += h;
		std::cout << "i: " << i << std::endl;
		double Rho_pert = ComputeAll(x_pert, myMesh, myCAND_OBJ_GRAD);//, i, h);
		//x_pert = x;
		//x_pert(i) -= h;
		//double Rho_pert_n = ComputeAll(x_pert, myMesh, myCAND_OBJ_GRAD);//, i, h);
		//df_dx(i) = (1/(2*h)) * (Rho_pert - Rho_pert_n);	
		df_dx(i) = (1/h) * (Rho_pert - Rho);	
		std::cout << "Rho: " << Rho << std::endl;
		std::cout << "Rho_pert: " << Rho_pert << std::endl;
	}

	std::cout << "df_dx_pert: ******************************" << df_dx << std::endl;
	std::cout << "df_dx_pert: ******************************\n" << df_dx.transpose() << std::endl;
	h = 0;
	el = 0;
	Rho = ComputeAll(x, myMesh, myCAND_OBJ_GRAD);//, el, h);//End Uncomment for FD */ 
	//std::cout << "dRho_dx:\n" << dRho_dx << std::endl;


	int opt_dim = 12; //Uncomment for Optimization
	nlopt::opt opt("LD_MMA", 12);
  	std::vector<double> lb(opt_dim);
	for (int i = 0; i < opt_dim; i++) {
		lb[i] = 0;
	}
  	//lb[0] = -HUGE_VAL; lb[1] = 0;
  	opt.set_lower_bounds(lb);
  	opt.set_min_objective(myvfunc, &myMesh);
  	//my_constraint_data data[2] = { {2,0}, {-1,1} };
  	//my_constraint_data data2[2] = { {2,0}, {-1,1} };
  	opt.add_inequality_constraint(myvconstraint, &myMesh, 1e-8);
  	//opt.add_inequality_constraint(myvconstraint, &data[1], 1e-8);
  	opt.set_xtol_rel(1e-4);

	// try setting an algorithm parameter: 
  	opt.set_param("inner_maxeval", 123);
  	if (opt.get_param("inner_maxeval", 1234) != 123 || opt.get_param("not a param", 1234) != 1234 ||
      	opt.num_params() != 1 || std::string(opt.nth_param(0)) != "inner_maxeval") {
    	std::cerr << "failed to retrieve nlopt parameter" << std::endl;
    	return EXIT_FAILURE;
  	}

	std::vector<double> x_init(opt_dim);
  	x_init[0] = 0.025; x_init[1] = 0.025; x_init[2] = 0.025; x_init[3] = 0.175; x_init[4] = 0.175; x_init[5] = 0.025;
	x_init[6] = 0.01; x_init[7] = 0.01; x_init[8] = 0.01; x_init[9] = 0.02; x_init[10] = 0.02; x_init[11] = 0.02;
	
  	double minf;

	try{
    	opt.optimize(x_init, minf);
    	std::cout << "found minimum at f(" << x_init[0] << "," << x_init[1] << "," << x_init[2] << "," << x_init[3] << "," << x_init[4] << "," << x_init[5] << "," << x_init[6] << "," << x_init[7]
				<< "," << x_init[8] << "," << x_init[9] << "," << x_init[10] << "," << x_init[11] << ") = "
              	<< std::setprecision(10) << minf <<std::endl;
    	return std::fabs(minf - 0.5443310474) < 1e-3 ? EXIT_SUCCESS : EXIT_FAILURE;
  	}
  	catch(std::exception &e) {
    	std::cerr << "nlopt failed: " << e.what() << std::endl;
    	return EXIT_FAILURE;
  	} //End Uncomment for Optimization
	
}
