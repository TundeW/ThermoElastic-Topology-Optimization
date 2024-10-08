#include <iostream>
#include <fstream>
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
#include "WriteOutputFunctions.h"
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
	myMesh->circDom = myCAND_OBJ_GRAD.circDom;
	myMesh->dcircDom_dx = myCAND_OBJ_GRAD.dcircDom_dx;
	myMesh->HF_VolFrac = myCAND_OBJ_GRAD.HF_VolFrac;
	myMesh->dHF_VolFracdx = myCAND_OBJ_GRAD.dHF_VolFracdx;
	myMesh->opt_iter = myMesh->opt_iter + 1;

	Eigen::MatrixXd dMaxStress_dx = myMesh->dMaxStress_dx;
	std::cout << "Design Variable: " << x_eigen << std::endl;
	std::cout << "Objectives: " << myMesh->MaxStress << std::endl;
	//std::cout << "Objectives gradient: \n" << std::endl;
    if (!grad.empty()) {
		for (int i = 0; i < x.size(); i++) {
			
			grad[i] = dMaxStress_dx(i);
			//std::cout << grad[i] << std::endl;
		} 	
    }
	
	//Write Objectives to file
	std::string filename = "Results/objective_data.csv";
	// Create an output filestream object for Objective data files
    std::ofstream objectiveFile(filename, std::ios::app);
	WriteObjective(objectiveFile, myMesh->opt_iter, myMesh->MaxStress);
	// Close the file
    objectiveFile.close();
    
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

void multi_constraint(unsigned m, double *result, unsigned n, const double* x, double* grad, void *data)
{	
	MESH *myMesh = reinterpret_cast<MESH*>(data);
	std::cout << "Constraint: " << myMesh->HeatTransferRate << std::endl;
	int conLen = myMesh->conLen;
	int opt_dim = myMesh->opt_dim;
	int circDomLen = myMesh->circDom.rows();

    for (int jj = 0; jj < conLen; ++jj) {
        if (jj < 1) {
			result[jj]=myMesh->HeatTransferRate;
		}
        else if (jj < circDomLen + 1) {
			result[jj]=myMesh->circDom(jj-1);
		}
		else if (jj < circDomLen + 2) {
			result[jj] = 0.07 - myMesh->HF_VolFrac; // Min HF Vol Frac
		} 
		else {
			result[jj] = myMesh->HF_VolFrac - 0.08; // Max HF Vol Frac
		}
    }
	std::cout << "Constraint2: " << myMesh->HeatTransferRate << std::endl;
	if (grad) {
	//if (!grad.empty()) {
		for (int j = 0; j < conLen; j++) {
			for (int i = 0; i < opt_dim; i++) {
				if (j < 1) {
					grad[i] = myMesh->dHeatTransferRate_dx(i);
					//std::cout << grad[i] << std::endl;
				}
				else if (j < circDomLen + 1){
					grad[j*opt_dim + i] = myMesh->dcircDom_dx(j-1,i);
				}
				else if (j < circDomLen + 2) {
					//std::cout << "InHere3: " << std::endl;
					//std::cout << "myMesh->dHF_VolFracdx.rows(): " << myMesh->dHF_VolFracdx.rows() << ", myMesh->dHF_VolFracdx.cols(): " << myMesh->dHF_VolFracdx.cols() << std::endl;
					grad[j*opt_dim + i] = - myMesh->dHF_VolFracdx(i); // Min HF Vol Frac Sensitivities
				}
				else {
					//std::cout << "InHere4: " << std::endl;
					grad[j*opt_dim + i] = myMesh->dHF_VolFracdx(i); // Max HF Vol Frac Sensitivities
				}
				
			}
		}
  	}

	//Write Constraints to file
	std::string filename = "Results/constraints_data.csv";
	// Create an output filestream object for constraint data files
    std::ofstream constraintFile(filename, std::ios::app);
	WriteConstraints(constraintFile, conLen, myMesh->opt_iter, result);
	// Close the file
    constraintFile.close();

    return;
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
	nelx = 32; nely = 24; nelz = 16;
	//nelx = 10; nely = 5; nelz = 5;
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
	int conLen = 1;
	int opt_dim = 12;
	int opt_iter = 0;

	MESH myMesh = {coordinates, nodemap, mesh_centers, nel, dim, conLen, opt_dim, opt_iter, L, H, B, nelx, nely, nelz, Area, he, Hmat};
	CAND_OBJ_GRAD myCAND_OBJ_GRAD;

	Eigen::VectorXd x(12);
    x << 0.025, 0.025, 0.025, 0.075, 0.175, 0.025, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02;
	//x << 0.025, 0.025, 0.01, 0.02;
	//x << 0.464914, 0.0890919, 0.025, 0.175, 1.075, 0, 0.670067, 0.01, 0.01, 1.757, 0.02, 0;
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


	opt_dim = 12; //Uncomment for Optimization
	nlopt::opt opt("LD_MMA", 12);
  	std::vector<double> lb(opt_dim);
	std::vector<double> ub(opt_dim);
	for (int i = 0; i < opt_dim; i++) {
		lb[i] = 0;
		if (i < opt_dim/2) { //Pipe Locations
			if (i) { //If i is even, x-coordinate
				ub[i] = L;
			}
			else { // If i is odd, y-coordinate
				ub[i] = H;
			}
		}
		else { //Pipe Radii
			ub[i] = H/2;
		}
	}
  	//lb[0] = -HUGE_VAL; lb[1] = 0;
  	opt.set_lower_bounds(lb);
	opt.set_upper_bounds(ub);
  	opt.set_min_objective(myvfunc, &myMesh);
  	//my_constraint_data data[2] = { {2,0}, {-1,1} };
  	//my_constraint_data data2[2] = { {2,0}, {-1,1} };
  	//opt.add_inequality_constraint(myvconstraint, &myMesh, 1e-8);
	conLen = opt_dim+3; //op_dim Constraints for Bounding the Circles, 1 Constraint for Physics (HT/Stress), 2 Constraints for HF Volume Fraction
	myMesh.conLen = conLen;
	myMesh.opt_dim = opt_dim;
	std::vector<double> tol_constraint(conLen,1e-8);
  	opt.add_inequality_mconstraint(multi_constraint, &myMesh, tol_constraint);
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
  	x_init[0] = 0.025; x_init[1] = 0.025; x_init[2] = 0.025; x_init[3] = 0.075; x_init[4] = 0.175; x_init[5] = 0.025;
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
