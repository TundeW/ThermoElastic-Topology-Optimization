#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseLU>
#include <iostream>
#include <chrono>
#include "ThermalFeaFunction.h"
#include <omp.h>

using namespace Eigen;

// Eigen::MatrixXd ThermalFea(const Eigen::MatrixXd Rho, const MESH myMesh, const Eigen::MatrixXd Hmat) {
Eigen::MatrixXd ThermalFea(const Eigen::MatrixXd Rho, const MESH myMesh, const Eigen::SparseMatrix<double> Hmat, FEA& myFea, Eigen::SparseMatrix<double>& Kg, GRAD_HELPER_VARS& myGrad_Helper_Vars) {
    std::cout << "Begin Thermal Analysis: " << std::endl;
	// Access struct parameters
    Eigen::MatrixXd coordinates = myMesh.coordinate;
    Eigen::MatrixXd nodemap = myMesh.nodemap;
    Eigen::MatrixXd mesh_centers = myMesh.mesh_centers;
	double nel = myMesh.nel;
	int dim = myMesh.dim;
	double nelx = myMesh.nelx;
	double nely = myMesh.nely;
	double nelz = myMesh.nelz;
	double he = myMesh.he;
	MatrixXd Rho_nodes = Hmat * Rho;
	//std::cout << "Rho_nodes:\n" << Rho_nodes << "\n\n";

	//Function to obtain the velocity distribution
	std::string velocity_filename = "velocity_data2.csv";
	//std::string velocity_filename = "velocity_data_tawk.csv";
	std::vector<std::vector<double>> velocity = read_csv(velocity_filename);

	//MatrixXd wt(1, pow(2, dim));
	//MatrixXd gPoints(pow(2, dim), dim);
	MatrixXd wt;
	MatrixXd gPoints;
	wt.resize(1, pow(2,dim));
	gPoints.resize(pow(2,dim), dim);

	if (dim == 2) {
		//wt.resize(1, 4);
		//gPoints.resize(4, 2);
		wt << 1, 1, 1, 1;
        gPoints << -1/sqrt(3), -1/sqrt(3),
                    1/sqrt(3), -1/sqrt(3),
                    1/sqrt(3),  1/sqrt(3),
                    -1/sqrt(3),  1/sqrt(3);
	}
	else if (dim == 3) {
		//wt.resize(1, 8);
		//gPoints.resize(8, 3);
		wt << 1, 1, 1, 1, 1, 1, 1, 1;
        gPoints << -1/sqrt(3), -1/sqrt(3), -1/sqrt(3),
                    1/sqrt(3), -1/sqrt(3), -1/sqrt(3),
                    1/sqrt(3),  1/sqrt(3), -1/sqrt(3),
                    -1/sqrt(3),  1/sqrt(3), -1/sqrt(3),
                    -1/sqrt(3), -1/sqrt(3),  1/sqrt(3),
                    1/sqrt(3), -1/sqrt(3),  1/sqrt(3),
                    1/sqrt(3),  1/sqrt(3),  1/sqrt(3),
            	    -1/sqrt(3),  1/sqrt(3),  1/sqrt(3);
	}

	int ndof = nodemap.maxCoeff() + 1;
	//std::cout << "ndof: " << ndof << std::endl;
	//MatrixXd Kg(ndof, ndof);
	//MatrixXd P(ndof, 1);
	//MatrixXd Kg;
	//Kg = MatrixXd::Zero(ndof, ndof); //Change 8 to 2**dim
	//***std::cout << "Creating Themal Stiffness Sparse Matrix" << std::endl;
	/*Eigen::SparseMatrix<double> Kg(ndof, ndof);
	Kg.reserve(VectorXi::Constant(ndof,27));

	std::vector<Eigen::MatrixXd> dKthdrho1(nel, Eigen::MatrixXd::Zero(pow(2,dim), pow(2,dim)));
	std::vector<Eigen::MatrixXd> dKthdrho2(nel, Eigen::MatrixXd::Zero(pow(2,dim), pow(2,dim)));
	std::vector<Eigen::MatrixXd> dPthdrho1(nel, Eigen::MatrixXd::Zero(pow(2,dim), 1));
	std::vector<Eigen::MatrixXd> dPthdrho2(nel, Eigen::MatrixXd::Zero(pow(2,dim), 1));*/

	//Kg = MatrixXd::Zero(8, 8); 
	MatrixXd P;
	P = MatrixXd::Zero(ndof, 1);
	MatrixXd Vel;
	Vel = MatrixXd::Zero(nel, dim);
	MatrixXd Kappa_Tau;
	Kappa_Tau = MatrixXd::Zero(nel, 3);
	//P = MatrixXd::Zero(8, 1);
	//MatrixXd T(ndof, 1);
	double p = 3;
	double kappa = 26.1; //
	double Q = 0; // Internal Heat generation
	double q = 0; // Surface flux
	double Tenv = 30; //Environment Temperature
	double hConv = 1; // Coefficient of Convection
	double rho_cp = 1;

	auto start_time = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < nel; i++) {

		MatrixXd enodes = nodemap.row(i);
		//std::cout << "enodes: " << enodes << std::endl;
		MatrixXd elem_coords(enodes.cols(), dim); // convert dim to int, dim=3
		for (int i = 0; i < enodes.cols(); i++) {
        	elem_coords.row(i) = coordinates.row(enodes(i));
			//std::cout << "elem_coords(" << i << "): " << elem_coords.row(i) << std::endl;
    	}

		
		double q = 0.01;
		double rho1Pen = Rho(i,0) * (1 + q)/(Rho(i,0) + q);
		double rho2Pen = Rho(i,1) * (1 + q)/(Rho(i,1) + q);
		double drho1Pen_drho1 = (q + pow(q, 2))/pow(Rho(i,0) + q, 2);
		double drho2Pen_drho2 = (q + pow(q, 2))/pow(Rho(i,1) + q, 2);
		//rho1Pen = 0.01 * (1 - Rho(i,0))/(Rho(i,0) + 0.01);
		//rho2Pen = 0.01 * (1 - Rho(i,1))/(Rho(i,1) + 0.01);
		double kf1 = 0.02; double kf2 = 0.02; double ks = 26.1;

		kappa = rho1Pen*(rho2Pen*kf2 + (1-rho2Pen)*kf1) + (1 - rho1Pen)*ks;
		double dkappa_drho1Pen = (rho2Pen*kf2 + (1-rho2Pen)*kf1) - ks;
		double dkappa_drho2Pen = rho1Pen*(kf2 - kf1);
		double dkappa_drho1 = dkappa_drho1Pen * drho1Pen_drho1;
		double dkappa_drho2 = dkappa_drho2Pen * drho2Pen_drho2;
		//kappa = rho1Pen*(rho2Pen*1.6611e-5 + (1-rho2Pen)*1.6611e-5) + (1 - rho1Pen)*0.0217e0;
		//kappa = rho1Pen*(rho2Pen*26.1 + (1-rho2Pen)*1) + (1 - rho1Pen)*26.1;
		//kappa = 0.04;
		rho_cp = 1.204 * 1000; //1.204Kg/m3 * 1000 J/KgK
		//rho_cp = 12.04;

		/*if (mesh_centers(i,0) > 0.46 & mesh_centers(i,0) < 0.54) {
			kappa = 20;
		}*/	

		double rho_min = 1e-3;
		double Rho1 = rho_min + (1-rho_min)*Rho(i,0);
    	double Rho2 = rho_min + (1-rho_min)*Rho(i,1);

		rho1Pen = Rho1 * (1 + 0.01)/(Rho1 + 0.01);
		rho2Pen = Rho2 * (1 + 0.01)/(Rho2 + 0.01);
		drho1Pen_drho1 = (q + pow(q, 2))/pow(Rho1 + q, 2) * (1-rho_min);
		drho2Pen_drho2 = (q + pow(q, 2))/pow(Rho2 + q, 2) * (1-rho_min);

		// Define Elemental Velocity Vector
		double v1 = 1;
		v1 = rho1Pen*(rho2Pen*1 + (1-rho2Pen)*0) + (1 - rho1Pen)*0; //Uncomment 
		double dv1_drho1Pen = (rho2Pen*1 + (1-rho2Pen)*0) - 0;
		double dv1_drho2Pen = rho1Pen*(1 - 0);
		double dv1_drho1 = dv1_drho1Pen * drho1Pen_drho1;
		double dv1_drho2 = dv1_drho2Pen * drho2Pen_drho2;

		double v2 = 0;
		double v3 = 0;
		v3 = rho1Pen*(rho2Pen*(0) + (1-rho2Pen)*(-1)) + (1 - rho1Pen)*0; //Uncomment 
		v3 = rho1Pen*(rho2Pen*(0) + (1-rho2Pen)*(1)) + (1 - rho1Pen)*0; 
		double dv3_drho1Pen = (rho2Pen*0 + (1-rho2Pen)*(1)) - 0;
		double dv3_drho2Pen = rho1Pen*(0 - 1);
		double dv3_drho1 = dv3_drho1Pen * drho1Pen_drho1;
		double dv3_drho2 = dv3_drho2Pen * drho2Pen_drho2;
		//v1 = velocity[i][0];
		//v2 = velocity[i][1];
		//v3 = velocity[i][2];
		MatrixXd V;
		V.resize(1, dim);
		V << v1, v2, v3;
		Vel.row(i) = V;
		double Vnorm = V.norm();

		MatrixXd dV_drho1;
		dV_drho1.resize(1, dim);
		dV_drho1 << dv1_drho1, v2, dv3_drho1;
		MatrixXd dV_drho2;
		dV_drho2.resize(1, dim);
		dV_drho2 << dv1_drho2, v2, dv3_drho2;
		//MatrixXd V_diag = V.asDiagonal();
		//Eigen::DiagonalMatrix<double, Eigen::Dynamic> V_diag(V);
		Eigen::DiagonalMatrix<double, Eigen::Dynamic> V_diag(V.transpose());

		// Define Elemental Material Properties
		double C_k = kappa/ks;
		//C_k = kappa;
		double U = 1;
		double L = 0.25; //L=he;
		double Pe = (rho_cp * U)/ks; 
		//double Pe = (rho_cp * Vnorm)/26.1; 
		//Pe = (rho_cp * U); 

		MatrixXd D_vect;
    	D_vect = MatrixXd::Constant(dim, 1, C_k);
		MatrixXd D;
    	D = D_vect.asDiagonal();

		MatrixXd dD_vect_drho1;
    	dD_vect_drho1 = MatrixXd::Constant(dim, 1, dkappa_drho1/ks);
		MatrixXd dD_drho1;
    	dD_drho1 = dD_vect_drho1.asDiagonal();

		MatrixXd dD_vect_drho2;
    	dD_vect_drho2 = MatrixXd::Constant(dim, 1, dkappa_drho2/ks);
		MatrixXd dD_drho2;
    	dD_drho2 = dD_vect_drho2.asDiagonal();

		//Stabilization Parameters
		double r = 2;
		double tau1 = 4*he/(Vnorm*Pe);
		//tau1 = 4*he/(U*Pe);
		double tau3 = he*he/(4*C_k);
		//tau3 = he*he/(4*26.1);
		double powmin = -1/r;
		double tau_SU = pow(1/pow(tau1,r) + 1/pow(tau3,r), powmin);
		Kappa_Tau(i,0) = kappa;
		Kappa_Tau(i,1) = tau_SU;

		MatrixXd dtau_SU_due = -pow(tau_SU, 3) * (Pe/(16*he*he)) * V; //Validate this later on!!!!!! 
		double dtau_SU_dtau1 = pow(tau1, -r-1) * pow(1/pow(tau1,r) + 1/pow(tau3,r), powmin-1);
		double dtau_SU_dtau3 = pow(tau3, -r-1) * pow(1/pow(tau1,r) + 1/pow(tau3,r), powmin-1);
		double V_V_transpose = (V*V.transpose())(0,0);
		double dtau1_drho1 = -tau1 * pow(V_V_transpose,-1) * (V * dV_drho1.transpose())(0,0);
		double dtau1_drho2 = -tau1 * pow(V_V_transpose,-1) * (V * dV_drho2.transpose())(0,0);
		double dtau3_drho1 = -(he*he/(4*C_k*C_k)) * dkappa_drho1/ks;
		double dtau3_drho2 = -(he*he/(4*C_k*C_k)) * dkappa_drho2/ks;
		double dtau_SU_drho1 = dtau_SU_dtau1*dtau1_drho1 + dtau_SU_dtau3*dtau3_drho1;
		double dtau_SU_drho2 = dtau_SU_dtau1*dtau1_drho2 + dtau_SU_dtau3*dtau3_drho2;


		MatrixXd k_el;
		k_el = MatrixXd::Zero(pow(2,dim), pow(2,dim)); 
		MatrixXd dk_el_drho1;
		dk_el_drho1 = MatrixXd::Zero(pow(2,dim), pow(2,dim)); 
		MatrixXd dk_el_drho2;
		dk_el_drho2 = MatrixXd::Zero(pow(2,dim), pow(2,dim)); 

		MatrixXd k_transport;
		k_transport = MatrixXd::Zero(pow(2,dim), pow(2,dim)); 
		MatrixXd dk_transport_drho1;
		dk_transport_drho1 = MatrixXd::Zero(pow(2,dim), pow(2,dim)); 
		MatrixXd dk_transport_drho2;
		dk_transport_drho2 = MatrixXd::Zero(pow(2,dim), pow(2,dim)); 

		MatrixXd k_supg;
		k_supg = MatrixXd::Zero(pow(2,dim), pow(2,dim)); //Change 8 to 2**dim
		MatrixXd dk_supg_drho1;
		dk_supg_drho1 = MatrixXd::Zero(pow(2,dim), pow(2,dim));
		MatrixXd dk_supg_drho2;
		dk_supg_drho2 = MatrixXd::Zero(pow(2,dim), pow(2,dim));
		//K_el = MatrixXd::Zero(8, 8); 
		MatrixXd p_el;
		p_el = MatrixXd::Zero(pow(2,dim), 1);
		//p_el = MatrixXd::Zero(8, 1);

		for (int a = 0; a < pow(2,dim); a++) {
			MatrixXd N;
			MatrixXd dN;
			shapeFcn(a, gPoints, myMesh, N, dN);
			//std::cout << "ipt_coords: \n" << ipt_coords << std::endl;

			MatrixXd J = elem_coords.transpose() * dN;
			double detJ = J.determinant();
			MatrixXd B_T = J.ldlt().solve(dN.transpose());
			MatrixXd B = B_T.transpose();

			k_el = k_el + B * D * B.transpose() * detJ * wt(a);
			dk_el_drho1 = dk_el_drho1 + B * dD_drho1 * B.transpose() * detJ * wt(a);
			dk_el_drho2 = dk_el_drho2 + B * dD_drho2 * B.transpose() * detJ * wt(a);

			k_transport = k_transport + N * V * B.transpose() * detJ * wt(a) * Pe;
			dk_transport_drho1 = dk_transport_drho1 + N * dV_drho1 * B.transpose() * detJ * wt(a) * Pe;
			dk_transport_drho2 = dk_transport_drho2 + N * dV_drho2 * B.transpose() * detJ * wt(a) * Pe;

			k_supg = k_supg + B * V.transpose() * V * B.transpose() * detJ * wt(a) * tau_SU * Pe * Pe;
			MatrixXd dk_supg_drho1_a = (B * ((V.transpose() * dV_drho1)+(dV_drho1.transpose() * V)) * B.transpose() * detJ * wt(a) * tau_SU * Pe * Pe) + (B * V.transpose() * V * B.transpose() * detJ * wt(a) * dtau_SU_drho1 * Pe * Pe);
			MatrixXd dk_supg_drho2_a = (B * ((V.transpose() * dV_drho2)+(dV_drho2.transpose() * V)) * B.transpose() * detJ * wt(a) * tau_SU * Pe * Pe) + (B * V.transpose() * V * B.transpose() * detJ * wt(a) * dtau_SU_drho2 * Pe * Pe);
			dk_supg_drho1 = dk_supg_drho1 + dk_supg_drho1_a;
			dk_supg_drho2 = dk_supg_drho2 + dk_supg_drho2_a;

			p_el = p_el + Q * N * detJ * wt(a);
		
		}

		// Matrix Assembly
		//MatrixXd enodes = nodemap.row(i);
		for (int ii = 0; ii < enodes.cols(); ++ii) {
			int enodes_i = enodes(ii);
    		for (int j = 0; j < enodes.cols(); ++j) {				
				int enodes_j = enodes(j);
				//Kg(enodes_i, enodes_j) = Kg(enodes_i, enodes_j) + k_el(ii, j);
				// Kg.coeffRef(enodes_i, enodes_j) = Kg.coeff(enodes_i, enodes_j) + k_el(ii, j);
				Kg.coeffRef(enodes_i, enodes_j) += (k_el(ii, j) + k_transport(ii, j) + k_supg(ii, j));

				//std::cout << "k_el(i, j): \n" << k_el(i, j) << std::endl;
				//std::cout << "Kg(enodes(i), enodes(j)): \n" << Kg(enodes_i, enodes_j) << std::endl;
    		}

    		// Load Vector Assembly
    		P(enodes_i) = P(enodes_i) + p_el(ii);
		}
		myGrad_Helper_Vars.dKthdrho1[i] += dk_el_drho1 + dk_transport_drho1 + dk_supg_drho1;
		myGrad_Helper_Vars.dKthdrho2[i] += dk_el_drho2 + dk_transport_drho2 + dk_supg_drho2;
		/*if (i==14) {
			std::cout << "dV_drho1: " << dV_drho1 << std::endl;
			std::cout << "dV_drho2: " << dV_drho2 << std::endl;
			std::cout << "dk_el_drho1: " << dk_el_drho1 << std::endl;
			std::cout << "dk_transport_drho1: " << dk_transport_drho1 << std::endl;
			std::cout << "dk_supg_drho1: " << dk_supg_drho1 << std::endl;
			std::cout << "dk_el_drho2: " << dk_el_drho2 << std::endl;
			std::cout << "dk_transport_drho2: " << dk_transport_drho2 << std::endl;
			std::cout << "dk_supg_drho2: " << dk_supg_drho2 << std::endl;
		}*/

		//myFea.K1 = 0;
		/*if (i==1) {
			MatrixXd enodes1 = nodemap.row(1);
			double K1_el_0_0 = k_supg(0,0);
			
			myFea.K1 = K1_el_0_0;
			std::cout << "K1_el_0_0: " << dtau1_drho2 << std::endl;
			std::cout << "dk_supg_drho2: " << dk_supg_drho2 << std::endl;
			std::cout << "myGrad_Helper_Vars.dKthdrho1[i]: " << myGrad_Helper_Vars.dKthdrho1[i] << std::endl;
			std::cout << "myGrad_Helper_Vars.dKthdrho2[i]: " << myGrad_Helper_Vars.dKthdrho2[i] << std::endl;
		}*/
	}

	/*std::cout << "myGrad_Helper_Vars.dKthdrho1[13]:\n" << myGrad_Helper_Vars.dKthdrho1[13] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho1[14]:\n" << myGrad_Helper_Vars.dKthdrho1[14] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho2[13]:\n" << myGrad_Helper_Vars.dKthdrho2[13] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho2[14]:\n" << myGrad_Helper_Vars.dKthdrho2[14] << std::endl;*/

	//std::cout << "dKthdrho1: \n" << dKthdrho1[0] << std::endl;

	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    //***std::cout << "Thermal Stiffness Assembly time: " << duration.count() << " seconds" << std::endl;
	//std::cout << "Kg: \n" << Kg << std::endl;

	//Boundary Conditions
	start_time = std::chrono::high_resolution_clock::now();
	//VectorXi fixedDofs = (coordinates.col(0).array() == 0.0).cast<int>();
	//Get Node data for fixed and free nodes
	MatrixXi fdof1;
	MatrixXi fdof2;

	for (int a = 0; a < coordinates.rows(); a++) {
		if (coordinates(a,0) == 0.0) {
			fdof1.conservativeResize(1, fdof1.cols() + 1);
    		fdof1.col(fdof1.cols() - 1) << a;
		}
		if (coordinates(a,0) == myMesh.L) {
		//if (coordinates(a,2) == 0.0) {
			fdof2.conservativeResize(1, fdof2.cols() + 1);
    		fdof2.col(fdof2.cols() - 1) << a;
		}
		/*if (coordinates(a,1) == 0.0 & coordinates(a,0) > 0.54) {
			fdof1.conservativeResize(1, fdof1.cols() + 1);
    		fdof1.col(fdof1.cols() - 1) << a;
		}
		if (coordinates(a,1) == myMesh.H & coordinates(a,0) < 0.46) {
		//if (coordinates(a,2) == 0.0) {
			fdof2.conservativeResize(1, fdof2.cols() + 1);
    		fdof2.col(fdof2.cols() - 1) << a;
		}*/
	}

	//Get Element data for Heat Convection and flux
	VectorXi x_y_front = VectorXi::LinSpaced(nelx*nely, 0, nelx*nely-1);
	VectorXi x_y_back = VectorXi::LinSpaced(nelx*nely, nelx*nely*(nelz-1), nelx*nely*nelz-1);
	VectorXi x__z_top = VectorXi::LinSpaced(nelx*nelz, 0, nelx*nely*nelz - nelz);
	VectorXi x_z_bottom = VectorXi::LinSpaced(nelx*nelz, nely - 1, nelx*nely*nelz - 1);

	VectorXi elemnum = VectorXi::LinSpaced(nely, 0, nely-1);
	VectorXi nodeidz = VectorXi::LinSpaced(nelz, 0, nelx*nely*(nelz-1));
	VectorXi elemdup = elemnum.replicate(nodeidz.rows(), nodeidz.cols());
	MatrixXi nodedup = nodeidz.replicate(elemnum.cols(), elemnum.rows());
	VectorXi nodedup_reshaped = nodedup.transpose().reshaped(elemdup.rows(), elemdup.cols());
	VectorXi y_z_left = elemdup + nodedup_reshaped;

	elemnum = VectorXi::LinSpaced(nely, (nelx-1)*nely, nelx*nely-1);
	elemdup = elemnum.replicate(nodeidz.rows(), nodeidz.cols());
	VectorXi y_z_right = elemdup + nodedup_reshaped;

	VectorXi x_y_front_nodes = VectorXi::LinSpaced(nelx*(nely+1), nely+1, (nelx+1)*(nely+1)-1);
	//fdof2 = x_y_front_nodes.transpose();
	//std::cout << "fdof2.cols(): " << fdof2.cols() << " fdof2:\n" << fdof2 << std::endl;
	//MatrixXi fdof1(1, 8);
	//MatrixXi fdof2(1, 8);
	//MatrixXi fixeddofs(1, 16);
	//MatrixXi fixeddofs(1, fdof1.cols() + fdof2.cols());
	MatrixXi fixeddofs(1, fdof1.cols()); //Uncomment 
	//MatrixXi freedofs(1, 16);
	//fdof1 << 0, 1, 2, 3, 16, 17, 18, 19;
	//fdof2 << 12, 13, 14, 15, 28, 29, 30, 31;
	//fixeddofs << fdof1, fdof2;
	fixeddofs << fdof1; //, fdof2; //Uncomment 
	//fixeddofs << 0, 1, 2, 3, 16, 17, 18, 19, 12, 13, 14, 15, 28, 29, 30, 31;
	//freedofs << 4, 5, 6, 7, 8, 9, 10, 11, 20, 21, 22, 23, 24, 25, 26, 27;
	Eigen::VectorXi sortedFixedDofs = fixeddofs.transpose();
    std::sort(sortedFixedDofs.data(), sortedFixedDofs.data() + sortedFixedDofs.size());


	VectorXi alldofs = VectorXi::LinSpaced(ndof, 0, ndof - 1);
	VectorXi freeDofs(alldofs.size());
	auto it = std::set_difference(alldofs.data(), alldofs.data() + alldofs.size(), 
                                  sortedFixedDofs.data(), sortedFixedDofs.data() + sortedFixedDofs.size(), 
                                  freeDofs.data());
    freeDofs.conservativeResize(std::distance(freeDofs.data(), it));
	MatrixXi freedofs = freeDofs.transpose();
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    //***std::cout << "Defining Boundary Conditions time: " << duration.count() << " seconds" << std::endl;
	
	// Convection 
	start_time = std::chrono::high_resolution_clock::now();
	// For y_z_right --> Change FluxEl and gPointsEdge for a different Edge 
	// VectorXi fluxEl = y_z_right;
	VectorXi fluxEl = x_y_back;
	//fluxEl = x_y_front;
	MatrixXd gPointsEdge(4, dim);
	/*gPointsEdge << 1, -1/sqrt(3), -1/sqrt(3),
               		1, 1/sqrt(3), -1/sqrt(3),
               		1,  -1/sqrt(3), 1/sqrt(3),
               		1,  1/sqrt(3), 1/sqrt(3);*/

	gPointsEdge << -1/sqrt(3), -1/sqrt(3), 1,
               		1/sqrt(3), -1/sqrt(3), 1,
               		1/sqrt(3), 1/sqrt(3), 1,
               	   -1/sqrt(3), 1/sqrt(3), 1;

	//std::cout << "fluxEl: " << fluxEl << std::endl;
	for (int a = 0; a < fluxEl.rows(); a++) {
		int el = fluxEl(a); // Element Number
		MatrixXd edge_nodes = nodemap.row(el);
		MatrixXd edge_elem_coords(edge_nodes.cols(), dim); // convert dim to int, dim=3
		for (int i = 0; i < edge_nodes.cols(); i++) {
        	edge_elem_coords.row(i) = coordinates.row(edge_nodes(i));
			//std::cout << "elem_coords(" << i << "): " << elem_coords.row(i) << std::endl;
    	}

		double w = (1 - Rho(el,1))*Rho(el,0);
		double dw_drho1 = 1 - Rho(el,1);
		double dw_drho2 = - Rho(el,0);
		//w=1;
		hConv = w/(1-w + 1e-6);
		double dhConv_dw = (1+1e-6)/pow(1-w+1e-6, 2);
		double dhConv_drho1 = dhConv_dw*dw_drho1;
		double dhConv_drho2 = dhConv_dw*dw_drho2;
		//std::cout << "w(" << a << ", " << el << "): " << w << ", hConv: " << hConv <<std::endl;	
		//hConv = 5;
		Tenv = 0; //1088.706; //0;
		double val = hConv * Tenv;
		double dval_drho1 = dhConv_drho1 * Tenv;
		double dval_drho2 = dhConv_drho2 * Tenv;

		MatrixXd wtEdge(1, 4);
		wtEdge << 1, 1, 1, 1; //Gauss Quadrature Weights

		MatrixXd kConv;
		kConv = MatrixXd::Zero(pow(2,dim), pow(2,dim)); //Change 8 to 2**dim
		MatrixXd dkConv_drho1;
		dkConv_drho1 = MatrixXd::Zero(pow(2,dim), pow(2,dim));
		MatrixXd dkConv_drho2;
		dkConv_drho2 = MatrixXd::Zero(pow(2,dim), pow(2,dim));

		MatrixXd pConv;
		pConv = MatrixXd::Zero(pow(2,dim), 1);
		MatrixXd dpConv_drho1;
		dpConv_drho1 = MatrixXd::Zero(pow(2,dim), 1);
		MatrixXd dpConv_drho2;
		dpConv_drho2 = MatrixXd::Zero(pow(2,dim), 1);

		MatrixXd pFlux;
		pFlux = MatrixXd::Zero(pow(2,dim), 1);

		for (int aa = 0; aa < pow(2,dim-1); aa++) {
			MatrixXd N;
			MatrixXd dN;
			shapeFcn(aa, gPointsEdge, myMesh, N, dN);
			//std::cout << "ipt_coords: \n" << ipt_coords << std::endl;
			MatrixXd J = edge_elem_coords.transpose() * dN;
			double dS = sqrt(pow(J(0,0),2) + pow(J(1,0),2)); //Check dS for 3D
			kConv = kConv + hConv * N * N.transpose() * dS * wtEdge(aa); //where hConv = 1
			dkConv_drho1 = dkConv_drho1 + dhConv_drho1 * N * N.transpose() * dS * wtEdge(aa);
			dkConv_drho2 = dkConv_drho2 + dhConv_drho2 * N * N.transpose() * dS * wtEdge(aa);

			pConv = pConv + val * N * dS * wtEdge(aa); 
			dpConv_drho1 = dpConv_drho1 + dval_drho1 * N * dS * wtEdge(aa); 
			dpConv_drho2 = dpConv_drho2 + dval_drho2 * N * dS * wtEdge(aa); 

			pFlux = pFlux + q * N * dS * wtEdge(aa); 
			/*if (el==14 & aa==0) {
				std::cout << "N: " << N << std::endl;
				std::cout << "dN: " << dN << std::endl;
				std::cout << "dS: " << dS << std::endl;
				std::cout << "J: " << J << std::endl;
				std::cout << "edge_elem_coords: " << edge_elem_coords << std::endl;
				//::cout << "Multiply: " << edge_elem_coords << std::endl;
			}*/
		}
		for (int ii = 0; ii < edge_nodes.cols(); ++ii) {
			int enodes_i = edge_nodes(ii);
    		for (int j = 0; j < edge_nodes.cols(); ++j) {
				int enodes_j = edge_nodes(j);
				//Kg(enodes_i, enodes_j) = Kg(enodes_i, enodes_j) + kConv(ii, j);
				Kg.coeffRef(enodes_i, enodes_j) +=  + kConv(ii, j); //Uncomment 
    		}
    		P(enodes_i) = P(enodes_i) + pConv(ii); // + pFlux(ii); //Uncomment 
		}

		/*if (el==14) {
			std::cout << "dhConv_drho1: " << dhConv_drho1 << std::endl;
			std::cout << "dhConv_drho2: " << dhConv_drho2 << std::endl;
			std::cout << "dkConv_drho1: " << dkConv_drho1 << std::endl;
			std::cout << "dkConv_drho2: " << dkConv_drho2 << std::endl;
		}*/

		myGrad_Helper_Vars.dKthdrho1[el] += dkConv_drho1;
		myGrad_Helper_Vars.dKthdrho2[el] += dkConv_drho2;
		myGrad_Helper_Vars.dPthdrho1[el] += dpConv_drho1;
		myGrad_Helper_Vars.dPthdrho2[el] += dpConv_drho2;

	}
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    //***std::cout << "Defining Convection Matrices time: " << duration.count() << " seconds" << std::endl;
	/*std::cout << "myGrad_Helper_Vars.dKthdrho1[12]:\n" << myGrad_Helper_Vars.dKthdrho1[12] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho1[13]:\n" << myGrad_Helper_Vars.dKthdrho1[13] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho1[14]:\n" << myGrad_Helper_Vars.dKthdrho1[14] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho2[12]:\n" << myGrad_Helper_Vars.dKthdrho2[12] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho2[13]:\n" << myGrad_Helper_Vars.dKthdrho2[13] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho2[14]:\n" << myGrad_Helper_Vars.dKthdrho2[14] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho2[0]:\n" << myGrad_Helper_Vars.dKthdrho2[0] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho1[2]:\n" << myGrad_Helper_Vars.dKthdrho1[2] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho2[2]:\n" << myGrad_Helper_Vars.dKthdrho2[2] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho1[20]:\n" << myGrad_Helper_Vars.dKthdrho1[20] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho2[20]:\n" << myGrad_Helper_Vars.dKthdrho2[20] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho1[23]:\n" << myGrad_Helper_Vars.dKthdrho1[23] << std::endl;
	std::cout << "myGrad_Helper_Vars.dKthdrho2[23]:\n" << myGrad_Helper_Vars.dKthdrho2[23] << std::endl;*/

	start_time = std::chrono::high_resolution_clock::now();
	VectorXi stiffEl = y_z_right;//x_y_front; //;
	/*gPointsEdge << -1/sqrt(3), -1/sqrt(3), -1,
            		1/sqrt(3), -1/sqrt(3), -1,
            		1/sqrt(3),  1/sqrt(3), -1, 
           			-1/sqrt(3),  1/sqrt(3), -1;
			
	gPointsEdge << 1, -1/sqrt(3), -1/sqrt(3),
               		1, 1/sqrt(3), -1/sqrt(3),
               		1,  -1/sqrt(3), 1/sqrt(3),
               		1,  1/sqrt(3), 1/sqrt(3);*/

	for (int a = 0; a < stiffEl.rows(); a++) {
		int el = fluxEl(a); // Element Number
		MatrixXd elem_nodes = nodemap.row(el);
		MatrixXd elem_coords(elem_nodes.cols(), dim); 
		for (int i = 0; i < elem_nodes.cols(); i++) {
        	elem_coords.row(i) = coordinates.row(elem_nodes(i));
			//std::cout << "elem_coords(" << i << "): " << elem_coords.row(i) << std::endl;
    	}
		MatrixXd kse;
		//kse = MatrixXd::Identity(pow(2,dim), pow(2,dim)); //Change 8 to 2**dim
		kse = MatrixXd::Zero(pow(2,dim), pow(2,dim)); 
		//double rho1 = Rho(el,1);
		//double rho2 = Rho(el,2);
		double stiff = 100; // * std::pow(rho1*(1-rho2),3)
		MatrixXd wtEdge(1, 4);
		wtEdge << 1, 1, 1, 1; //Gauss Quadrature Weights
		//kse = stiff * MatrixXd::Identity(pow(2,dim), pow(2,dim));

		for (int aa = 0; aa < pow(2,dim-1); aa++) {
			MatrixXd N;
			MatrixXd dN;
			shapeFcn(aa, gPointsEdge, myMesh, N, dN);
			//std::cout << "ipt_coords: \n" << ipt_coords << std::endl;
			MatrixXd J = elem_coords.transpose() * dN;
			double dS = sqrt(pow(J(1,1),2) + pow(J(2,1),2));
			kse = kse + stiff * N * N.transpose() * dS * wtEdge(aa); //where hConv = 1
			
		}
		//std::cout << "Kse:\n" << kse << std::endl;
		for (int ii = 0; ii < elem_nodes.cols(); ++ii) {
			int enodes_i = elem_nodes(ii);
    		for (int j = 0; j < elem_nodes.cols(); ++j) {
				int enodes_j = elem_nodes(j);
				// Kg(enodes_i, enodes_j) = Kg(enodes_i, enodes_j) + kse(ii, j);
    		}
		}
	}
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    //***std::cout << "Defining Elemental stiffness penalization time: " << duration.count() << " seconds" << std::endl;
	
	start_time = std::chrono::high_resolution_clock::now();
	//VectorXi x_y_front_nodes = VectorXi::LinSpaced(nelx*(nely+1), nely+1, (nelx+1)*(nely+1)-1);
	for (int a = 0; a < x_y_front_nodes.rows(); a++) {
		int node_i = x_y_front_nodes(a);
		double rho1Pen = Rho_nodes(node_i,0) * (1 + 0.01)/(Rho_nodes(node_i,0) + 0.01);
		double rho2Pen = Rho_nodes(node_i,1) * (1 + 0.01)/(Rho_nodes(node_i,1) + 0.01);
		//rho1Pen = 0.01 * (1 - Rho_nodes(node_i,0))/(Rho_nodes(node_i,0) + 0.01);
		//rho2Pen = 0.01 * (1 + Rho_nodes(node_i,1))/(Rho_nodes(node_i,1) + 0.01);

		double stiff = rho1Pen*(rho2Pen*0 + (1-rho2Pen)*1000) + (1 - rho1Pen)*0;

		// stiff = 100;
		// Kg(node_i, node_i) = Kg(node_i, node_i) + stiff;
		// Kg.coeffRef(node_i, node_i) = Kg.coeff(node_i, node_i) + stiff;
		// Kg.coeffRef(node_i, node_i) += stiff;
	}
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    //***std::cout << "Defining Nodal stiffness penalization time: " << duration.count() << " seconds" << std::endl;
	
	start_time = std::chrono::high_resolution_clock::now();
	//std::cout << "Rho size: " << Rho.rows() << "x" << Rho.cols() << std::endl;
	// Create Temperature Vector
	MatrixXd T;
	//T = MatrixXd::Zero(ndof, 1);
	T = MatrixXd::Zero(ndof, 1);

	for (int a = 0; a < fdof1.cols(); a++) {
		//T(fdof2(a)) = 0; //Uncomment 
		T(fdof1(a)) = 100; //533.15; //100;	 //Uncomment 
		//T(fdof2(a)) = 1;  
		//T(fdof1(a)) = 0;	 
	}
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    //***std::cout << "Creating T matrix: " << duration.count() << " seconds" << std::endl;

	start_time = std::chrono::high_resolution_clock::now();
	//Solve the thermal Problem Matrix System
	MatrixXd T_fixeddofs(fixeddofs.cols(), 1);
	MatrixXd P_freedofs(freedofs.cols(), 1);
	// MatrixXd Kg_freefree(freedofs.cols(), freedofs.cols());
	Eigen::SparseMatrix<double> Kg_freefree(freedofs.cols(), freedofs.cols());
	Kg_freefree.reserve(VectorXi::Constant(freedofs.cols(),27));
	// MatrixXd Kg_freefixed(freedofs.cols(), fixeddofs.cols());
	Eigen::SparseMatrix<double> Kg_freefixed(freedofs.cols(), fixeddofs.cols());
	Kg_freefixed.reserve(VectorXi::Constant(fixeddofs.cols(),27));

	//std::cout << "fixed.cols & free.cols()" << fixeddofs.cols() << " & " << freedofs.cols() << std::endl;
	//Eigen::SparseMatrix<double> Kg(ndof, ndof);
	//Kg.reserve(VectorXi::Constant(ndof,27));

	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    //***std::cout << "Initializing Matrix for free and fixed dofs indexing: " << duration.count() << " seconds" << std::endl;

	//int numThreads = Eigen::nbThreads();
    //std::cout << "Number of threads: " << numThreads << std::endl;
	//Eigen::setNbThreads(6);
	//std::cout << "Number of threads: " << numThreads << std::endl;
	start_time = std::chrono::high_resolution_clock::now();
	//#pragma omp parallel for
	for (int b = 0; b < freedofs.cols(); b++) {
		P_freedofs(b) = P(freedofs(b));
		//Kg_freefixed.col(b) = Kg.coeff(freedofs(b), fixeddofs).transpose(); //Potentially use alternative to loops to improve performance
    	//Kg_freefree.col(b) = Kg.coeff(freedofs(b), freedofs).transpose();
		for (int a = 0; a < fixeddofs.cols(); a++) {
			T_fixeddofs(a) = T(fixeddofs(a));
			//#pragma omp critical
			if (Kg.coeff(freedofs(b), fixeddofs(a)) != 0.0) 
			{
				Kg_freefixed.insert(b,a) = Kg.coeff(freedofs(b), fixeddofs(a));
			}
		}
		for (int c = 0; c < freedofs.cols(); c++) {

			//#pragma omp critical
			if (Kg.coeff(freedofs(b), freedofs(c)) != 0.0) 
			{
				Kg_freefree.insert(b,c) = Kg.coeff(freedofs(b), freedofs(c));
			}
		}
	}
	myGrad_Helper_Vars.Kth_freefree = Kg_freefree;
	myGrad_Helper_Vars.th_freedofs = freedofs;
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    //***std::cout << "Inseting values into the Kg and T fixed and free dofs: " << duration.count() << " seconds" << std::endl;

	start_time = std::chrono::high_resolution_clock::now();
	MatrixXd Pg = P_freedofs - Kg_freefixed * T_fixeddofs;
	// MatrixXd T_freedofs = Kg_freefree.colPivHouseholderQr().solve(Pg);
	// MatrixXd T_freedofs = Kg_freefree.ldlt().solve(Pg);
	

	

	//SparseLU<SparseMatrix<double> > solver;
	BiCGSTAB<SparseMatrix<double> > solver;

	 // Create the Incomplete LU preconditioner
    /*Eigen::IncompleteLUT<double> preconditioner;
    preconditioner.compute(Kg_freefree);  // Compute the ILU factorization
    // Set the preconditioner for the solver
    solver.preconditioner() = preconditioner;*/

	//Solve the system
	solver.compute(Kg_freefree);
	MatrixXd T_freedofs = solver.solve(Pg);
	Eigen::setNbThreads(Eigen::nbThreads());  // Reset the number of threads to the default

	if (solver.info() != Eigen::Success) {
        //std::out << "SparseLU solver failed to factorize the matrix!" << std::endl;
		std::cout << "BiCGSTAB solver failed to converge!" << std::endl;
        // return 1;
    }
	for (int c = 0; c < freedofs.cols(); c++) {
		T(freedofs(c)) = T_freedofs(c);
	}	
	/*std::cout << "T(" << freedofs(freedofs.cols()-1) << ")" << T(freedofs(freedofs.cols()-1)) << std::endl;
	std::cout << "T(" << freedofs(freedofs.cols()-2) << ")" << T(freedofs(freedofs.cols()-2)) << std::endl;
	std::cout << "T(" << freedofs(freedofs.cols()-3) << ")" << T(freedofs(freedofs.cols()-3)) << std::endl;
	std::cout << "T(" << freedofs(freedofs.cols()-4) << ")" << T(freedofs(freedofs.cols()-4)) << std::endl;
	std::cout << "T(" << freedofs(freedofs.cols()-5) << ")" << T(freedofs(freedofs.cols()-5)) << std::endl;*/
	
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    //***std::cout << "Solving Linear System time: " << duration.count() << " seconds" << std::endl;
	myFea.Vel = Vel;
	myFea.Kappa_Tau = Kappa_Tau;

	/**********************************For Finite Difference Validation of Kg************************************************/
	MatrixXd enodes1 = nodemap.row(17);
	int enodes_i = enodes1(7);
	int enodes_j = enodes1(7);
	double K1_el_0_0 = Kg.coeff(enodes_i, enodes_j);
	double P1_el_0_0 = P(enodes_i);
	myFea.K1 = K1_el_0_0;
	//std::cout << "K1_el_0_0: " << K1_el_0_0 << std::endl;
	/**********************************************************************************/
	
	
	return T;
}