#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <chrono>
#include "ElasticFeaFunction.h"

using namespace Eigen;

Eigen::MatrixXd ElasticFea(const Eigen::MatrixXd Rho, const MESH myMesh, const Eigen::MatrixXd Hmat, const Eigen::MatrixXd T, FEA& myFea, Eigen::SparseMatrix<double>& Kg, GRAD_HELPER_VARS& myGrad_Helper_Vars) {
	std::cout << "Begin Elastic Analysis: " << std::endl;
	// Access struct parameters
    Eigen::MatrixXd coordinates = myMesh.coordinate;
    Eigen::MatrixXd nodemap = myMesh.nodemap;
    Eigen::MatrixXd mesh_centers = myMesh.mesh_centers;
	double nel = myMesh.nel;
	int dim = myMesh.dim;
	double nelx = myMesh.nelx;
	double nely = myMesh.nely;
	double nelz = myMesh.nelz;
	MatrixXd Rho_nodes = Hmat * Rho;
	//std::cout << "Rho_nodes:\n" << Rho_nodes << "\n\n";

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

	int ndof = dim * (nodemap.maxCoeff() + 1);
	//MatrixXd Kg;
	//Kg = MatrixXd::Zero(ndof, ndof); //Change 8 to 2**dim
	//***std::cout << "Creating Elastic Stiffness Sparse Matrix" << std::endl;
	//Eigen::SparseMatrix<double> Kg(ndof, ndof);
	//Kg.reserve(VectorXi::Constant(ndof,81));

	MatrixXd P;
	P = MatrixXd::Zero(ndof, 1);

	//Defining derivative Info that needs to be accessed outside the function
	/*Eigen::SparseMatrix<double> dFdT(ndof, ndof/dim);
	dFdT.reserve(VectorXi::Constant(81,ndof/dim));
	std::vector<Eigen::MatrixXd> dKeldrho1(nel, Eigen::MatrixXd::Zero(pow(2,dim), pow(2,dim)));
	std::vector<Eigen::MatrixXd> dKeldrho2(nel, Eigen::MatrixXd::Zero(pow(2,dim), pow(2,dim)));
	std::vector<Eigen::MatrixXd> dPeldrho1(nel, Eigen::MatrixXd::Zero(pow(2,dim), 1));
	std::vector<Eigen::MatrixXd> dPeldrho2(nel, Eigen::MatrixXd::Zero(pow(2,dim), 1));
	std::vector<Eigen::MatrixXd> dstrain_thermaldrho1(nel, Eigen::MatrixXd::Zero(6, 1));
	std::vector<Eigen::MatrixXd> dstrain_thermaldrho2(nel, Eigen::MatrixXd::Zero(6, 1));
	std::vector<Eigen::MatrixXd> dstrain_thermaldt(nel, Eigen::MatrixXd::Zero(6, pow(2,dim)));*/

	double p = 3; //Interpolation Penalization
	double E = 1.75e+11; //175 GPA: dynamic modulus of elasticity of haynes 282 @ 800 degree C
	double nu = 0.355; //Poisson ratio at 800 degree C
	double alpha=14.9e-6; //14.9 µm/m.°C: mean coefficient of thermal expansion of haynes 282 @800 degrees C
	double Tenv = 30; //Environment Temperature

	MatrixXd D;
	MatrixXd strain_thermal;
	MatrixXd A;

	if (dim == 2) {
		D.resize(3, 3);
		A.resize(3, 1);
		//gPoints.resize(4, 2);
		strain_thermal = Eigen::MatrixXd::Zero(3, nel);
        D << 1, nu, 0,
             nu, 1, 0,
             0,  0, (1-nu)/2;
		A << 1, 1, 0;
                    
	}
	else if (dim == 3) {
		D.resize(6, 6);
		A.resize(6, 1);
		strain_thermal = Eigen::MatrixXd::Zero(6, nel);
        D << 1-nu, nu, nu, 0, 0, 0,
             nu, 1-nu, nu, 0, 0, 0,
             nu, nu, 1-nu, 0, 0, 0,
             0,  0,  0, 0.5-nu, 0, 0,
             0,  0,  0, 0, 0.5-nu, 0,
             0,  0,  0, 0, 0, 0.5-nu;
		A << 1, 1, 1, 0, 0, 0;
	}

	auto start_time = std::chrono::high_resolution_clock::now();
	MatrixXd DD = 1/((1+nu)*(1-2*nu)) * D;
	myFea.D_Elastic = DD;
	for (int i = 0; i < nel; i++) {

		MatrixXd enodes = nodemap.row(i);
		MatrixXd elem_coords(enodes.cols(), dim); // convert dim to int, dim=3
		for (int ii = 0; ii < enodes.cols(); ii++) {
        	elem_coords.row(ii) = coordinates.row(enodes(ii));
    	}

		//Elastic Modulus
		double Emin = E*1e-9;
		double E_rho = Emin + pow((1-Rho(i,0)),3)*(E - Emin);
		double dE_rho_drho1 = -3*pow((1-Rho(i,0)),2)*(E - Emin);
		double dE_rho_drho2 = 0;

		myFea.Kappa_Tau(i,2) = E_rho;
		MatrixXd D_rho = E_rho*DD;
		MatrixXd dD_rho_drho1 = dE_rho_drho1*DD;
		MatrixXd dD_rho_drho2 = dE_rho_drho2*DD;

		//Thermal Strain
		MatrixXd cenPoint = Eigen::MatrixXd::Zero(1, 3);
		MatrixXd N;
		MatrixXd dN;
		shapeFcn(0, cenPoint, myMesh, N, dN);
		MatrixXd T_el(enodes.cols(), 1); // convert dim to int, dim=3


		for (int ii = 0; ii < enodes.cols(); ii++) {
			int enodes_i = enodes(ii);
        	T_el(ii) = T(enodes_i);
    	}
		MatrixXd Te_vect = N.transpose() * T_el;;
		double Te = Te_vect(0,0);
		double dTe = Te - Tenv;

		/*//Strain Thermal
		double rho = 1e-3 + (1 - 1e-3)*(1-Rho(i,0));
		strain_thermal.col(i) = pow(rho,3) * alpha * dTe * A;
		myGrad_Helper_Vars.dstrain_thermaldrho1[i] = - (1 - 1e-3) * 3 * pow(rho,2) * alpha * dTe * A;
		myGrad_Helper_Vars.dstrain_thermaldrho2[i] = 0 * alpha * dTe * A;
		myGrad_Helper_Vars.dstrain_thermaldt[i] = pow(rho,3) * alpha * A * N.transpose();*/

		// Alternative Strain Thermal
		strain_thermal.col(i) = pow((1-Rho(i,0)),3) * alpha * dTe * A;
		myGrad_Helper_Vars.dstrain_thermaldrho1[i] = - 3 * pow((1-Rho(i,0)),2) * alpha * dTe * A;
		myGrad_Helper_Vars.dstrain_thermaldrho2[i] = 0 * alpha * dTe * A;
		myGrad_Helper_Vars.dstrain_thermaldt[i] = pow((1-Rho(i,0)),3) * alpha * A * N.transpose();
		

		/*if (dim == 2) {
			strain_thermal.col(i) << alpha*dTe, alpha*dTe, 0;
                    
		}
		else if (dim == 3) {
			strain_thermal.col(i) << alpha*dTe, alpha*dTe, alpha*dTe, 0, 0, 0;
			
		}*/

		MatrixXd k_el;
		MatrixXd f_thermal_el;
		MatrixXd df_thermal_eldt;
		MatrixXd df_thermal_eldstrain_thermal;
		MatrixXd dk_el_drho1;
		MatrixXd dk_el_drho2;
		MatrixXd delf_thermal_el_delrho1;
		MatrixXd delf_thermal_el_delrho2;

		elem_stiff(D_rho, elem_coords, gPoints, strain_thermal.col(i), wt, myMesh.Area, myMesh, k_el, f_thermal_el, df_thermal_eldt, df_thermal_eldstrain_thermal, myGrad_Helper_Vars.dstrain_thermaldt[i], dk_el_drho1, dk_el_drho2, dD_rho_drho1, dD_rho_drho2, delf_thermal_el_delrho1, delf_thermal_el_delrho2);
		//elem_stiff(D****, elem_coords, gPoints, strain_thermal*******, wt, **********A, myMesh, k_el, f_thermal_el, df_thermal_eldt, df_thermal_eldstrain_thermal, dstrain_thermaldt***, dk_el_drho1, dk_el_drho2, dD_rho_drho1, dD_rho_drho2);

		MatrixXd edof = nodes_to_dofs(enodes, myMesh); 

		for (int ii = 0; ii < edof.cols(); ++ii) {
			int edof_i = edof(ii);
    		for (int j = 0; j < edof.cols(); ++j) {
				int edof_j = edof(j);
				//Kg(edof_i, edof_j) = Kg(edof_i, edof_j) + k_el(ii, j);
				Kg.coeffRef(edof_i, edof_j) += k_el(ii, j);
			}
    		P(edof_i) = P(edof_i) + f_thermal_el(ii);
			for (int j = 0; j < enodes.cols(); ++j) {
				int enodes_j = enodes(j);
				myGrad_Helper_Vars.dFdT.coeffRef(edof_i, enodes_j) += df_thermal_eldt(ii,j);
			}
		}
		myGrad_Helper_Vars.dKeldrho1[i] = dk_el_drho1;
		myGrad_Helper_Vars.dKeldrho2[i] = dk_el_drho2;
		myGrad_Helper_Vars.dPeldrho1[i] = df_thermal_eldstrain_thermal*myGrad_Helper_Vars.dstrain_thermaldrho1[i] + delf_thermal_el_delrho1;
		myGrad_Helper_Vars.dPeldrho2[i] = df_thermal_eldstrain_thermal*myGrad_Helper_Vars.dstrain_thermaldrho2[i] + delf_thermal_el_delrho2;

		/*if (i == 0 || i == 2 || i == 8 || i == 14) {
			std::cout << "strain_thermal.col(i)(" << i << "): " << strain_thermal.col(i) << std::endl;
			std::cout << "dTe: " << dTe << std::endl;
			std::cout << "df_thermal_eldstrain_thermal(" << i << "): " << df_thermal_eldstrain_thermal << std::endl;
			std::cout << "myGrad_Helper_Vars.dstrain_thermaldrho1(" << i << "): " << myGrad_Helper_Vars.dstrain_thermaldrho1[i] << std::endl;
		}*/

	}
	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    //***std::cout << "Elastic Stiffness Matrix time: " << duration.count() << " seconds" << std::endl;
    
	start_time = std::chrono::high_resolution_clock::now();
	myFea.strain_thermal = strain_thermal;
	MatrixXi fnodes1;
	MatrixXi fnodes2;

	for (int a = 0; a < coordinates.rows(); a++) {
		if (coordinates(a,2) == 0.0) {
			fnodes1.conservativeResize(1, fnodes1.cols() + 1);
    		fnodes1.col(fnodes1.cols() - 1) << a;
		}
		if (coordinates(a,2) == myMesh.B) {
			fnodes2.conservativeResize(1, fnodes2.cols() + 1);
    		fnodes2.col(fnodes2.cols() - 1) << a;
		}
	}

	MatrixXi fixednodes(1, fnodes1.cols() + fnodes2.cols());
	fixednodes << fnodes1, fnodes2;

	// ************************************************************ //
	// Validation Test
	/*for (int a = 0; a < coordinates.rows(); a++) {
		if (coordinates(a,0) == 0.0) {
			fnodes1.conservativeResize(1, fnodes1.cols() + 1);
    		fnodes1.col(fnodes1.cols() - 1) << a;
		}
		if (coordinates(a,0) == myMesh.L) {
			fnodes2.conservativeResize(1, fnodes2.cols() + 1);
    		fnodes2.col(fnodes2.cols() - 1) << a;
			//P(3*a + 1) = -1;
		}
		//if (coordinates(a,0) == myMesh.L & coordinates(a,1) == myMesh.H ) {
		//	fnodes2.conservativeResize(1, fnodes2.cols() + 1);
    	//	fnodes2.col(fnodes2.cols() - 1) << a;
			//P(3*a + 1) = -1;
		//}
	}

	//std::cout << "The Load Nodes: " << fnodes2 << std::endl;

	//MatrixXi fixednodes(1, fnodes1.cols());
	//fixednodes << fnodes1;
	MatrixXi fixednodes(1, fnodes1.cols() + fnodes2.cols());
	fixednodes << fnodes1, fnodes2;*/
	// *************************************************************** //

	MatrixXd fixednodes_double = fixednodes.cast<double>();
	MatrixXd fixeddofs_double = nodes_to_dofs(fixednodes_double, myMesh);
	MatrixXi fixeddofs = fixeddofs_double.cast<int>();

	Eigen::VectorXi sortedFixedDofs = fixeddofs.transpose();
    std::sort(sortedFixedDofs.data(), sortedFixedDofs.data() + sortedFixedDofs.size());


	VectorXi alldofs = VectorXi::LinSpaced(ndof, 0, ndof - 1);
	VectorXi freeDofs(alldofs.size());
	auto it = std::set_difference(alldofs.data(), alldofs.data() + alldofs.size(), 
                                  sortedFixedDofs.data(), sortedFixedDofs.data() + sortedFixedDofs.size(), 
                                  freeDofs.data());
    freeDofs.conservativeResize(std::distance(freeDofs.data(), it));
	MatrixXi freedofs = freeDofs.transpose();
	//std::cout << "fixeddofs: " << fixeddofs.rows() << "x" << fixeddofs.cols() << std::endl;
	//std::cout << "freedofs: " << freedofs.rows() << "x" << freedofs.cols() << std::endl;

	MatrixXd U = MatrixXd::Zero(ndof, 1);
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    //***std::cout << "Elastic Boundary Conditions time: " << duration.count() << " seconds" << std::endl;
    

	/*for (int a = 0; a < fixeddofs.cols(); a++) {
		U(fixeddofs(a)) = 0;
		//T(fdof1(a)) = 100;
		
	}*/
	
	start_time = std::chrono::high_resolution_clock::now();
	//Solve the Elastic Problem Matrix System
	MatrixXd U_fixeddofs(fixeddofs.cols(), 1);
	MatrixXd P_freedofs(freedofs.cols(), 1);
	//MatrixXd Kg_freefree(freedofs.cols(), freedofs.cols());
	//MatrixXd Kg_freefixed(freedofs.cols(), fixeddofs.cols());
	Eigen::SparseMatrix<double> Kg_freefree(freedofs.cols(), freedofs.cols());
	Kg_freefree.reserve(VectorXi::Constant(freedofs.cols(),81));
	Eigen::SparseMatrix<double> Kg_freefixed(freedofs.cols(), fixeddofs.cols());
	Kg_freefixed.reserve(VectorXi::Constant(fixeddofs.cols(),81));
	
	/* // For Dense Matrices
	for (int b = 0; b < freedofs.cols(); b++) {
		P_freedofs(b) = P(freedofs(b));
		for (int a = 0; a < fixeddofs.cols(); a++) {
			U_fixeddofs(a) = U(fixeddofs(a));
			Kg_freefixed(b,a) = Kg(freedofs(b), fixeddofs(a));
		}
		for (int c = 0; c < freedofs.cols(); c++) {
			Kg_freefree(b,c) = Kg(freedofs(b), freedofs(c));
		}
	}*/

	// For Sparse Matrices
	for (int b = 0; b < freedofs.cols(); b++) {
		P_freedofs(b) = P(freedofs(b));
		for (int a = 0; a < fixeddofs.cols(); a++) {
			U_fixeddofs(a) = U(fixeddofs(a));
			if (Kg.coeff(freedofs(b), fixeddofs(a)) != 0.0) 
			{
				Kg_freefixed.insert(b,a) = Kg.coeff(freedofs(b), fixeddofs(a));
			}
		}
		for (int c = 0; c < freedofs.cols(); c++) {
			if (Kg.coeff(freedofs(b), freedofs(c)) != 0.0) 
			{
				Kg_freefree.insert(b,c) = Kg.coeff(freedofs(b), freedofs(c));
			}
		}
	}
	myGrad_Helper_Vars.Kel_freefree = Kg_freefree;
	myGrad_Helper_Vars.el_freedofs = freedofs;
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    //***std::cout << "Inseting values into the Kg and U fixed and free dofs:" << duration.count() << " seconds" << std::endl;
	
	start_time = std::chrono::high_resolution_clock::now();
	//MatrixXd Pg = P_freedofs - Kg_freefixed * U_fixeddofs;
	// MatrixXd U_freedofs = Kg_freefree.colPivHouseholderQr().solve(P_freedofs);
	//MatrixXd U_freedofs = Kg_freefree.ldlt().solve(P_freedofs);

	//Sparse
	BiCGSTAB<SparseMatrix<double> > solver;
	solver.compute(Kg_freefree);
	MatrixXd U_freedofs = solver.solve(P_freedofs);
	Eigen::setNbThreads(Eigen::nbThreads());  // Reset the number of threads to the default

	for (int c = 0; c < freedofs.cols(); c++) {
		U(freedofs(c)) = U_freedofs(c);
	}
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    //***std::cout << "Solving Elastic Linear System time: " << duration.count() << " seconds" << std::endl;
    

	return U;
	
}