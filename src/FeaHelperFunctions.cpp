#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <fstream>
#include "FeaHelperFunctions.h"

using namespace Eigen;

void shapeFcn(const int i, const Eigen::MatrixXd gPoints, const MESH myMesh, MatrixXd& N, MatrixXd& dN) {
    // Access struct parameters
	int dim = myMesh.dim;
	// Define N and dN outside the if statement, using dim
	N.resize(pow(2,dim), 1);
    dN.resize(pow(2,dim), dim);
	MatrixXd ipt_coords;
	if (dim == 2) {

	}
	else if (dim == 3) {
		ipt_coords = gPoints.row(i);
		//MatrixXd NN(8, 1);
    	//MatrixXd dNN(8, 3);
		//N.resize(8, 1);
    	//dN.resize(8, 3);
		N << (1-ipt_coords(0))*(1-ipt_coords(1))*(1-ipt_coords(2))/8,
             (1+ipt_coords(0))*(1-ipt_coords(1))*(1-ipt_coords(2))/8,
             (1+ipt_coords(0))*(1+ipt_coords(1))*(1-ipt_coords(2))/8,
             (1-ipt_coords(0))*(1+ipt_coords(1))*(1-ipt_coords(2))/8,
             (1-ipt_coords(0))*(1-ipt_coords(1))*(1+ipt_coords(2))/8,
             (1+ipt_coords(0))*(1-ipt_coords(1))*(1+ipt_coords(2))/8,
             (1+ipt_coords(0))*(1+ipt_coords(1))*(1+ipt_coords(2))/8,
             (1-ipt_coords(0))*(1+ipt_coords(1))*(1+ipt_coords(2))/8;

		dN << -(1-ipt_coords(1))*(1-ipt_coords(2))/8, -(1-ipt_coords(0))*(1-ipt_coords(2))/8, -(1-ipt_coords(0))*(1-ipt_coords(1))/8,
               (1-ipt_coords(1))*(1-ipt_coords(2))/8, -(1+ipt_coords(0))*(1-ipt_coords(2))/8, -(1+ipt_coords(0))*(1-ipt_coords(1))/8,
               (1+ipt_coords(1))*(1-ipt_coords(2))/8,  (1+ipt_coords(0))*(1-ipt_coords(2))/8, -(1+ipt_coords(0))*(1+ipt_coords(1))/8,
              -(1+ipt_coords(1))*(1-ipt_coords(2))/8,  (1-ipt_coords(0))*(1-ipt_coords(2))/8, -(1-ipt_coords(0))*(1+ipt_coords(1))/8,
              -(1-ipt_coords(1))*(1+ipt_coords(2))/8, -(1-ipt_coords(0))*(1+ipt_coords(2))/8,  (1-ipt_coords(0))*(1-ipt_coords(1))/8,
               (1-ipt_coords(1))*(1+ipt_coords(2))/8, -(1+ipt_coords(0))*(1+ipt_coords(2))/8,  (1+ipt_coords(0))*(1-ipt_coords(1))/8,
               (1+ipt_coords(1))*(1+ipt_coords(2))/8,  (1+ipt_coords(0))*(1+ipt_coords(2))/8,  (1+ipt_coords(0))*(1+ipt_coords(1))/8,
              -(1+ipt_coords(1))*(1+ipt_coords(2))/8,  (1-ipt_coords(0))*(1+ipt_coords(2))/8,  (1-ipt_coords(0))*(1+ipt_coords(1))/8;

		//double coeff = 1/8;
		//N = coeff * NN;
		//dN = coeff * dNN;
		//std::cout << "ipt_coords: \n" << ipt_coords << std::endl;
		//std::cout << "N: \n" << N << std::endl;
		//std::cout << "dN: \n" << dN << std::endl;
	}
	//std::cout << "ipt_coords: \n" << ipt_coords << std::endl;
	//std::cout << "N: \n" << N << std::endl;
	//std::cout << "dN: \n" << dN << std::endl;
}

void elem_stiff(const Eigen::MatrixXd D, const Eigen::MatrixXd elem_coords, const Eigen::MatrixXd gPoints, const Eigen::MatrixXd strain_thermal, const Eigen::MatrixXd wt, const double A, const MESH myMesh, MatrixXd& k_el, MatrixXd& f_thermal_el, MatrixXd& df_thermal_eldt, MatrixXd& df_thermal_eldstrain_thermal, const Eigen::MatrixXd dstrain_thermaldt, MatrixXd& dk_el_drho1, MatrixXd& dk_el_drho2, const Eigen::MatrixXd dD_rho_drho1, const Eigen::MatrixXd dD_rho_drho2, MatrixXd& delf_thermal_el_delrho1, MatrixXd& delf_thermal_el_delrho2) {
	int dim = myMesh.dim;
	k_el = MatrixXd::Zero(dim*pow(2,dim), dim*pow(2,dim)); 
	dk_el_drho1 = MatrixXd::Zero(dim*pow(2,dim), dim*pow(2,dim)); 
	dk_el_drho2 = MatrixXd::Zero(dim*pow(2,dim), dim*pow(2,dim)); 

	f_thermal_el = MatrixXd::Zero(dim*pow(2,dim), 1); 

	df_thermal_eldt = MatrixXd::Zero(dim*pow(2,dim), pow(2,dim)); 
	df_thermal_eldstrain_thermal = MatrixXd::Zero(dim*pow(2,dim), 6); 
	delf_thermal_el_delrho1 = MatrixXd::Zero(dim*pow(2,dim), 1); 
	delf_thermal_el_delrho2 = MatrixXd::Zero(dim*pow(2,dim), 1); 

	for (int a = 0; a < pow(2,dim); a++) {
		MatrixXd N;
		MatrixXd dN;
		shapeFcn(a, gPoints, myMesh, N, dN);

		MatrixXd J = elem_coords.transpose() * dN;
		double detJ = J.determinant();
		MatrixXd dN_xy = J.ldlt().solve(dN.transpose());
		
		MatrixXd B;
		if (dim == 2) {
			B.resize(3, dim*pow(2,dim));
			for (int j = 0; j < pow(2,dim); j++) {
				B(0,2*j) = dN_xy(0,j);
                B(1,2*j+1) = dN_xy(1,j);
                B(2,2*j) = dN_xy(1,j);
                B(2,2*j+1) = dN_xy(0,j);
			}
		}
		else if (dim == 3) {
			B.resize(6, dim*pow(2,dim));
			B = MatrixXd::Zero(6, dim*pow(2,dim));
			for (int j = 0; j < pow(2,dim); j++) {
				B(0,3*j) = dN_xy(0,j);
                B(1,3*j+1) = dN_xy(1,j);
                B(2,3*j+2) = dN_xy(2,j);
                B(3,3*j+1) = dN_xy(2,j);
                B(3,3*j+2) = dN_xy(1,j);
                B(4,3*j) = dN_xy(2,j);
                B(4,3*j+2) = dN_xy(0,j);
                B(5,3*j)   = dN_xy(1,j);
                B(5,3*j+1) = dN_xy(0,j);
			}
		}

		k_el += B.transpose() * D * B * detJ * wt(a);
		dk_el_drho1 += B.transpose() * dD_rho_drho1 * B * detJ * wt(a);
		dk_el_drho2 += B.transpose() * dD_rho_drho2 * B * detJ * wt(a);

		f_thermal_el += A * detJ * B.transpose() * D * strain_thermal;
		df_thermal_eldt += A * detJ * B.transpose() * D * dstrain_thermaldt;
		df_thermal_eldstrain_thermal += A * detJ * B.transpose() * D;
		delf_thermal_el_delrho1 += A * detJ * B.transpose() * dD_rho_drho1 * strain_thermal;
		delf_thermal_el_delrho2 += A * detJ * B.transpose() * dD_rho_drho2 * strain_thermal;

		/*if ( a == 0 ) {
			std::cout << "dD_rho_drho1: " << dD_rho_drho1 << std::endl;
			std::cout << "dN_xy: " << dN_xy << std::endl;
			std::cout << "B: " << B << std::endl;
			std::cout << "detJ: " << detJ << std::endl;
			std::cout << "wt(a): " << wt(a) << std::endl;
		}*/
	}
}

Eigen::MatrixXd nodes_to_dofs(const Eigen::MatrixXd enodes, const MESH myMesh) {
	int n = enodes.cols();
	MatrixXd el_dofs = MatrixXd::Zero(1, myMesh.dim*n);

	if (myMesh.dim == 2) {
		for (int i = 0; i < n; i++) {
			el_dofs(2*i) = 2*enodes(i);
			el_dofs(2*i + 1) = 2*enodes(i) + 1;
		}
	}
	else if (myMesh.dim == 3) {
		for (int i = 0; i < n; i++) {
			el_dofs(3*i) = 3*enodes(i);
            el_dofs(3*i + 1) = 3*enodes(i) + 1;
            el_dofs(3*i + 2) = 3*enodes(i) + 2;
		}
	}
	return el_dofs;
}

double ComputeStress(const Eigen::MatrixXd Rho, const MESH myMesh, const Eigen::MatrixXd Hmat, FEA& myFea, GRAD_HELPER_VARS& myGrad_Helper_Vars) {
	std::cout << "Computing VonMises Stress at the Element Centroids" << std::endl;
	// Access struct parameters
    Eigen::MatrixXd nodemap = myMesh.nodemap;
	double nel = myMesh.nel;
	int dim = myMesh.dim;
	//MatrixXd Rho_nodes = Hmat * Rho;
	int ndof = dim * (nodemap.maxCoeff() + 1);
	MatrixXd strain_thermal = myFea.strain_thermal;
	MatrixXd D = myFea.D_Elastic;
	MatrixXd U = myFea.U;
	MatrixXd T = myFea.T;
	double pnorm = 10;

	MatrixXd strain_true;
	MatrixXd stress_true;
	MatrixXd vonMises;
	MatrixXd M;
	MatrixXd dfdu = Eigen::MatrixXd::Zero(ndof, 1);
	MatrixXd AdjEl = Eigen::MatrixXd::Zero(ndof, 1);
	MatrixXd Adjth = Eigen::MatrixXd::Zero(ndof/dim, 1);
	MatrixXd dfdT = Eigen::MatrixXd::Zero(ndof/dim, 1);
	MatrixXd delfdelrho1 = Eigen::MatrixXd::Zero(nel, 1);
	MatrixXd delfdelrho2 = Eigen::MatrixXd::Zero(nel, 1);

	if (dim == 2) {
		strain_true = Eigen::MatrixXd::Zero(3, nel);
		stress_true = Eigen::MatrixXd::Zero(3, nel);
		vonMises = Eigen::MatrixXd::Zero(nel, 1);
		M.resize(3, 3);
        M << 1, -0.5, 0,
             -0.5, 1, 0,
             0, 0, 3;
	}
	else if (dim == 3) {
		strain_true = Eigen::MatrixXd::Zero(6, nel);
		stress_true = Eigen::MatrixXd::Zero(6, nel);
		vonMises = Eigen::MatrixXd::Zero(nel, 1);
		M.resize(6, 6);
        M << 1, -0.5, -0.5, 0, 0, 0,
             -0.5, 1, -0.5, 0, 0, 0,
             -0.5, -0.5, 1, 0, 0, 0,
			 0,  0,  0, 3, 0, 0,
			 0,  0,  0, 0, 3, 0,
			 0,  0,  0, 0, 0, 3;
	}

	for (int i = 0; i < nel; i++) {
		MatrixXd enodes = nodemap.row(i);
		MatrixXd elem_coords(enodes.cols(), dim); // convert dim to int, dim=3
		for (int ii = 0; ii < enodes.cols(); ii++) {
        	elem_coords.row(ii) = myMesh.coordinate.row(enodes(ii));
    	}
		MatrixXd edof = nodes_to_dofs(enodes, myMesh); 

		MatrixXd U_elem(edof.cols(), 1);
		for (int ii = 0; ii < edof.cols(); ii++) {
			int edof_i = edof(ii);
        	U_elem(ii) = U(edof_i);
    	}

		//Elastic Modulus
		double E = 1.75e+11;
		double Emin = E*1e-9;
		double E_rho = Emin + pow((1-Rho(i,0)),3)*(E - Emin);
		MatrixXd D_rho = E_rho*D;
		double dE_rhodrho1 = -3 * pow((1-Rho(i,0)),2)*(E - Emin);
		MatrixXd dD_rhodrho1 = dE_rhodrho1 * D;


		MatrixXd cenPoint = Eigen::MatrixXd::Zero(1, 3); //Different for 2D
		MatrixXd strain_total = Eigen::MatrixXd::Zero(6, 1); //Different for 2D
		MatrixXd N;
		MatrixXd dN;
		shapeFcn(0, cenPoint, myMesh, N, dN);
		MatrixXd J = elem_coords.transpose() * dN;
		MatrixXd dN_xy = J.ldlt().solve(dN.transpose());
		
		MatrixXd B;
		if (dim == 2) {
			B.resize(3, dim*pow(2,dim));
			for (int j = 0; j < pow(2,dim); j++) {
				B(0,2*j) = dN_xy(0,j);
                B(1,2*j+1) = dN_xy(1,j);
                B(2,2*j) = dN_xy(1,j);
                B(2,2*j+1) = dN_xy(0,j);
			}
		}
		else if (dim == 3) {
			B.resize(6, dim*pow(2,dim));
			B = MatrixXd::Zero(6, dim*pow(2,dim));
			for (int j = 0; j < pow(2,dim); j++) {
				B(0,3*j) = dN_xy(0,j);
                B(1,3*j+1) = dN_xy(1,j);
                B(2,3*j+2) = dN_xy(2,j);
                B(3,3*j+1) = dN_xy(2,j);
                B(3,3*j+2) = dN_xy(1,j);
                B(4,3*j) = dN_xy(2,j);
                B(4,3*j+2) = dN_xy(0,j);
                B(5,3*j)   = dN_xy(1,j);
                B(5,3*j+1) = dN_xy(0,j);
			}
		}
		strain_total = B * U_elem;
		strain_true.col(i) = strain_total - strain_thermal.col(i);
		//strain_true.col(i) = strain_total;
		stress_true.col(i) = D_rho * strain_true.col(i);
		double rho = 1e-3 + (1 - 1e-3)*(1-Rho(i,0));
		vonMises(i) = pow(rho,3) * strain_true.col(i).transpose() * D_rho * M * D_rho * strain_true.col(i);

		//MatrixXd stress = Eigen::MatrixXd::Zero(3, nel);
		
		MatrixXd DfDu = (pnorm * pow(vonMises(i),pnorm-1) * 2 * pow(rho,3) * B.transpose() * D_rho * M * D_rho * strain_true.col(i));
		for (int ii = 0; ii < edof.cols(); ++ii) {
			int edof_i = edof(ii);
    		dfdu(edof_i) = dfdu(edof_i) + DfDu(ii);
		}

		MatrixXd DfDT = pnorm * pow(vonMises(i),pnorm-1) * (-2) * pow(rho,3) * myGrad_Helper_Vars.dstrain_thermaldt[i].transpose() * D_rho * M * D_rho * strain_true.col(i);
		for (int ii = 0; ii < enodes.cols(); ++ii) {
			int enodes_i = enodes(ii);
    		dfdT(enodes_i) = dfdT(enodes_i) + DfDT(ii);
		}
		//std::cout << "Rho2(" << i << "): " << Rho(i,1) << std::endl;
		//std::cout << "rho: " << rho << std::endl;

		MatrixXd alpha = -(1 - 1e-3) * 3 * pow(rho,2) * strain_true.col(i).transpose() * D_rho * M * D_rho * strain_true.col(i);
		MatrixXd beta = pow(rho,3) * (-2) * myGrad_Helper_Vars.dstrain_thermaldrho1[i].transpose() * D_rho * M * D_rho * strain_true.col(i);
		MatrixXd gamma = pow(rho,3) * 2 * strain_true.col(i).transpose() * dD_rhodrho1 * M * D_rho * strain_true.col(i);
		delfdelrho1(i) = pnorm * pow(vonMises(i),pnorm-1) * ( alpha(0,0) + beta(0,0) + gamma(0,0) );

		if (i == 0 || i == 2 || i == 8 || i == 14) {
			/*std::cout << "vonMises(i)(" << i << "): " << vonMises(i) << std::endl;
			std::cout << "Rho(i,0): " << Rho(i,0) << std::endl;
			std::cout << "rho: " << rho << std::endl;
			std::cout << "D_rho: " << D_rho << std::endl;
			std::cout << "M: " << M << std::endl;
			std::cout << "dE_rhodrho1: " << dE_rhodrho1 << std::endl;
			std::cout << "dD_rhodrho1(" << i << "): " << dD_rhodrho1 << std::endl;
			std::cout << "strain_true.col(i)(" << i << "): " << strain_true.col(i) << std::endl;
			std::cout << "B: " << B << std::endl;
			std::cout << "U_elem: " << U_elem << std::endl;
		//}
		//if (i==2) {
		std::cout << "DfDu: " << DfDu << std::endl;
		std::cout << "DfDT: " << DfDT << std::endl;
		std::cout << "myGrad_Helper_Vars.dstrain_thermaldt[i].transpose(): " << myGrad_Helper_Vars.dstrain_thermaldt[i].transpose() << std::endl;
		
		std::cout << "alpha_incomp: " << -(1 - 1e-3) * 3 * pow(rho,2) << std::endl;
		std::cout << "alpha: " << alpha << std::endl;
		std::cout << "beta: " << beta << std::endl;
		std::cout << "gamma: " << gamma << std::endl;
		std::cout << "delfdelrho1(i): " << delfdelrho1(i) << std::endl;
		std::cout << "D_rho: " << D_rho << std::endl;
		std::cout << "M: " << M << std::endl;
		std::cout << "strain_true.col(i): " << strain_true.col(i) << std::endl;*/

		/*std::cout << "Rho12[" << i << "]: " << Rho(i,0) << ", " << Rho(i,1) << std::endl;
		std::cout << "myGrad_Helper_Vars.dstrain_thermaldrho1[" << i << "]: " << myGrad_Helper_Vars.dstrain_thermaldrho1[i] << std::endl;
		std::cout << "myGrad_Helper_Vars.dstrain_thermaldrho2[" << i << "]: " << myGrad_Helper_Vars.dstrain_thermaldrho2[i] << std::endl;
		std::cout << "myGrad_Helper_Vars.dstrain_thermaldt[" << i << "]: " << myGrad_Helper_Vars.dstrain_thermaldt[i] << std::endl;
		std::cout << "rho: " << rho << std::endl;
		std::cout << "D_rho: " << D_rho << std::endl;
		std::cout << "M: " << M << std::endl;
		std::cout << "strain_true.col(i): " << strain_true.col(i) << std::endl;
		std::cout << "Product: " << myGrad_Helper_Vars.dstrain_thermaldrho1[i].transpose() * D_rho * M * D_rho * strain_true.col(i) << std::endl;
		MatrixXd Product = D_rho * M * D_rho * strain_true.col(i);
		std::cout << "Product: " << Product << std::endl;
		std::cout << "Product2: " << myGrad_Helper_Vars.dstrain_thermaldrho1[i].transpose() * Product << std::endl;
		std::cout << "beta: " << beta << std::endl;*/
		}
		
		
	}

	//std::cout << "dfdT: " << dfdT.transpose() << std::endl;
	//std::cout << "dfdu: " << dfdu.transpose() << std::endl;

	//objective
	double f = vonMises.lpNorm<10>();
	MatrixXd vonMises_pow = vonMises.array().pow(pnorm);
	//std::cout << "vonMises_pow: " << vonMises_pow << ", " << vonMises_pow.sum() << ", " << pow(vonMises_pow.sum(), (1/pnorm)-1) << ", " << pow(vonMises_pow.sum(), (1/pnorm)-1) * (1/pnorm) << std::endl;
	//std::cout << "(1/pnorm): " << 1/pnorm << std::endl;
	double phi = (1/pnorm) * pow(vonMises_pow.sum(), (1/pnorm)-1);
	//std::cout << "delfdelrho1: " << delfdelrho1 << std::endl;
	//std::cout << "phi: " << phi << std::endl;
	delfdelrho1 = phi * delfdelrho1;
	//std::cout << "delfdelrho1: " << delfdelrho1 << std::endl;
	MatrixXd df_du = phi * dfdu.transpose();
	MatrixXd df_dT = phi * dfdT.transpose();

	//std::cout << "vonMises: " << vonMises << std::endl;
	std::cout << "f: " << f << std::endl;
	//std::cout << "vonMises_pow: " << vonMises_pow << std::endl;
	//std::cout << "phi: " << phi << std::endl;

	/*std::cout << "df_dT: " << df_dT << std::endl;
	std::cout << "df_du: " << df_du << std::endl;*/

	//MatrixXd dRf_elduf = myGrad_Helper_Vars.Kel_freefree;
	Eigen::SparseMatrix<double> dRf_elduf = myGrad_Helper_Vars.Kel_freefree;
	//MatrixXd dRf_thdTf = myGrad_Helper_Vars.Kth_freefree;
	Eigen::SparseMatrix<double> dRf_thdTf = myGrad_Helper_Vars.Kth_freefree;
	MatrixXi th_freedofs = myGrad_Helper_Vars.th_freedofs;
	MatrixXi el_freedofs = myGrad_Helper_Vars.el_freedofs;

	MatrixXd dfduf = Eigen::MatrixXd::Zero(el_freedofs.cols(), 1);
	for (int c = 0; c < el_freedofs.cols(); c++) {
		dfduf(c) = df_du(el_freedofs(c));
	}

	Eigen::SparseMatrix<double> dReldTf(ndof, th_freedofs.cols());
	dReldTf.reserve(VectorXi::Constant(th_freedofs.cols(),81)); //*****Switched the rows and cols
	MatrixXd df_dTf = Eigen::MatrixXd::Zero(1, th_freedofs.cols());

	for (int c = 0; c < th_freedofs.cols(); c++) {
		df_dTf(c) = df_dT(th_freedofs(c));
		for (int b = 0; b < ndof; b++) {
			if (myGrad_Helper_Vars.dFdT.coeff(b, th_freedofs(c)) != 0.0) 
			{
				dReldTf.insert(b,c) = - myGrad_Helper_Vars.dFdT.coeff(b, th_freedofs(c));
			}
		}
	}


	//MatrixXd AdjEl_f = - dRf_elduf.ldlt().solve(dfduf);
	BiCGSTAB<SparseMatrix<double> > solver;
	solver.compute(dRf_elduf);
	MatrixXd AdjEl_f = solver.solve(-dfduf);
	for (int c = 0; c < el_freedofs.cols(); c++) {
		AdjEl(el_freedofs(c)) = AdjEl_f(c);
	}

	//std::cout << "AdjEl_f:\n" << AdjEl_f << std::endl;

	MatrixXd RHS = (AdjEl.transpose()*dReldTf + df_dTf);
	//MatrixXd Adjth_f = dRf_thdTf.ldlt().solve(RHS.transpose());
	BiCGSTAB<SparseMatrix<double> > solver2;
	solver2.compute(dRf_thdTf);
	MatrixXd Adjth_f = solver2.solve(-RHS.transpose());
	/*std::cout << "Size of AdjEl: (" << AdjEl.rows() << ", " << AdjEl.cols() << ")" << std::endl;
	std::cout << "Size of dReldTf: (" << dReldTf.rows() << ", " << dReldTf.cols() << ")" << std::endl;
	std::cout << "Size of df_dTf: (" << df_dTf.rows() << ", " << df_dTf.cols() << ")" << std::endl;*/
	for (int c = 0; c < th_freedofs.cols(); c++) {
		Adjth(th_freedofs(c)) = Adjth_f(c);
	}
	/*std::cout << "Adjth_f:\n" << Adjth_f << std::endl;
	std::cout << "dRf_thdTf(0,0):\n" << dRf_thdTf.coeff(0,0) << std::endl;
	std::cout << "df_dTf:\n" << df_dTf << std::endl;
	std::cout << "dReldTf(6,2):\n" << dReldTf.coeff(6,2) << std::endl;*/

	myGrad_Helper_Vars.dfdrho1.resize(1, nel);
	myGrad_Helper_Vars.dfdrho2.resize(1, nel);

	for (int i = 0; i < nel; ++i) {
		MatrixXd enodes = nodemap.row(i);
		MatrixXd edof = nodes_to_dofs(enodes, myMesh); 

		MatrixXd U_elem(edof.cols(), 1);
		MatrixXd AdjEl_elem(edof.cols(), 1);
		for (int ii = 0; ii < edof.cols(); ii++) {
			int edof_i = edof(ii);
        	U_elem(ii) = U(edof_i);
			AdjEl_elem(ii) = AdjEl(edof_i);
    	}

		MatrixXd T_el(enodes.cols(), 1); 
		MatrixXd Adjth_el(enodes.cols(), 1); 
		for (int ii = 0; ii < enodes.cols(); ii++) {
			int enodes_i = enodes(ii);
        	T_el(ii) = T(enodes_i);
			Adjth_el(ii) = Adjth(enodes_i);
    	}

		myGrad_Helper_Vars.dfdrho1(i) = delfdelrho1(i) + (AdjEl_elem.transpose()*(myGrad_Helper_Vars.dKeldrho1[i]*U_elem - myGrad_Helper_Vars.dPeldrho1[i]) + Adjth_el.transpose()*(myGrad_Helper_Vars.dKthdrho1[i]*T_el - myGrad_Helper_Vars.dPthdrho1[i]))(0,0);
		myGrad_Helper_Vars.dfdrho2(i) = delfdelrho2(i) + (AdjEl_elem.transpose()*(myGrad_Helper_Vars.dKeldrho2[i]*U_elem - myGrad_Helper_Vars.dPeldrho2[i]) + Adjth_el.transpose()*(myGrad_Helper_Vars.dKthdrho2[i]*T_el - myGrad_Helper_Vars.dPthdrho2[i]))(0,0);

		/*if (i==1) {
			std::cout << "delfdelrho2(i): " << delfdelrho2(i) << std::endl;
			std::cout << "dReldrho(i): " << AdjEl_elem.transpose()*(myGrad_Helper_Vars.dKeldrho2[i]*U_elem - myGrad_Helper_Vars.dPeldrho2[i]) << std::endl;
			std::cout << "dRthdrho(i): " << Adjth_el.transpose()*(myGrad_Helper_Vars.dKthdrho2[i]*T_el - myGrad_Helper_Vars.dPthdrho2[i]) << std::endl;
			std::cout << "myGrad_Helper_Vars.dfdrho2(i): " << myGrad_Helper_Vars.dfdrho2(i) << std::endl;
			std::cout << "myGrad_Helper_Vars.dKthdrho2[i]: " << myGrad_Helper_Vars.dKthdrho2[i] << std::endl;
			std::cout << "myGrad_Helper_Vars.dPthdrho2[i]: " << myGrad_Helper_Vars.dPthdrho2[i] << std::endl;
		}*/
		
	}

	//std::cout << "myGrad_Helper_Vars.dfdrho1: " << delfdelrho1 << std::endl;
	//std::cout << "myGrad_Helper_Vars.dfdrho2: " << delfdelrho2 << std::endl;

	myFea.stress_true = stress_true;
	myFea.vonMises = vonMises;
	return f;
}

double ComputeHeatTransferRate(const Eigen::MatrixXd Rho, const MESH myMesh, FEA& myFea, GRAD_HELPER_VARS& myGrad_Helper_Vars) {
	std::cout << "Computing Heat Transfer Rate" << std::endl;
	// Access struct parameters
    Eigen::MatrixXd coordinates = myMesh.coordinate;
    Eigen::MatrixXd nodemap = myMesh.nodemap;
    Eigen::MatrixXd mesh_centers = myMesh.mesh_centers;
	double nel = myMesh.nel;
	int dim = myMesh.dim;
	double nelx = myMesh.nelx;
	double nely = myMesh.nely;
	double nelz = myMesh.nelz;
	MatrixXi th_freedofs = myGrad_Helper_Vars.th_freedofs;
	
	//MatrixXd Rho_nodes = Hmat * Rho;
	int ndof = 1 * (nodemap.maxCoeff() + 1);
	MatrixXd strain_thermal = myFea.strain_thermal;

	MatrixXd T = myFea.T;
	
	Eigen::SparseMatrix<double> K_gFlux_inlet(ndof, ndof);
	K_gFlux_inlet.reserve(VectorXi::Constant(ndof,27));
	std::vector<Eigen::MatrixXd> dK_gFlux_inlet_drho1(nel, Eigen::MatrixXd::Zero(pow(2,dim), pow(2,dim)));
	std::vector<Eigen::MatrixXd> dK_gFlux_inlet_drho2(nel, Eigen::MatrixXd::Zero(pow(2,dim), pow(2,dim)));
	MatrixXd L_inlet = Eigen::MatrixXd::Zero(ndof, 1);

	Eigen::SparseMatrix<double> K_gFlux_outlet(ndof, ndof);
	K_gFlux_outlet.reserve(VectorXi::Constant(ndof,27));
	std::vector<Eigen::MatrixXd> dK_gFlux_outlet_drho1(nel, Eigen::MatrixXd::Zero(pow(2,dim), pow(2,dim)));
	std::vector<Eigen::MatrixXd> dK_gFlux_outlet_drho2(nel, Eigen::MatrixXd::Zero(pow(2,dim), pow(2,dim)));
	MatrixXd L_outlet = Eigen::MatrixXd::Zero(ndof, 1);

	VectorXi x_y_front = VectorXi::LinSpaced(nelx*nely, 0, nelx*nely-1);
	VectorXi x_y_back = VectorXi::LinSpaced(nelx*nely, nelx*nely*(nelz-1), nelx*nely*nelz-1);

	VectorXi inletEl = x_y_back;
	for (int i = 0; i < inletEl.rows(); i++) {
		int el = inletEl(i);
		MatrixXd enodes = nodemap.row(el);
		MatrixXd elem_coords(enodes.cols(), dim); // convert dim to int, dim=3
		for (int ii = 0; ii < enodes.cols(); ii++) {
        	elem_coords.row(ii) = coordinates.row(enodes(ii));
			int enodes_i = enodes(ii);
			L_inlet(enodes_i) = 1;
    	}

		double w = (1 - Rho(el,1))*Rho(el,0);
		double dw_drho1 = 1 - Rho(el,1);
		double dw_drho2 = - Rho(el,0);

		MatrixXd wtEdge(1, 4);
		wtEdge << 1, 1, 1, 1; //Gauss Quadrature Weights

		MatrixXd kConv;
		kConv = MatrixXd::Zero(pow(2,dim), pow(2,dim)); //Change 8 to 2**dim
		MatrixXd dkConv_drho1;
		dkConv_drho1 = MatrixXd::Zero(pow(2,dim), pow(2,dim));
		MatrixXd dkConv_drho2;
		dkConv_drho2 = MatrixXd::Zero(pow(2,dim), pow(2,dim));

		MatrixXd gPointsEdge(4, dim);

		gPointsEdge << -1/sqrt(3), -1/sqrt(3), 1,
               			1/sqrt(3), -1/sqrt(3), 1,
               			1/sqrt(3), 1/sqrt(3), 1,
               	   		-1/sqrt(3), 1/sqrt(3), 1;

		for (int aa = 0; aa < pow(2,dim-1); aa++) {
			MatrixXd N;
			MatrixXd dN;
			shapeFcn(aa, gPointsEdge, myMesh, N, dN);
			//std::cout << "ipt_coords: \n" << ipt_coords << std::endl;
			MatrixXd J = elem_coords.transpose() * dN;
			double dS = sqrt(pow(J(0,0),2) + pow(J(1,0),2)); //Check dS for 3D
			kConv = kConv + w * N * N.transpose() * dS * wtEdge(aa); //where hConv = 1
			dkConv_drho1 = dkConv_drho1 + dw_drho1 * N * N.transpose() * dS * wtEdge(aa);
			dkConv_drho2 = dkConv_drho2 + dw_drho2 * N * N.transpose() * dS * wtEdge(aa);

		}

		for (int ii = 0; ii < enodes.cols(); ++ii) {
			int enodes_i = enodes(ii);
    		for (int j = 0; j < enodes.cols(); ++j) {
				int enodes_j = enodes(j);
				//Kg(enodes_i, enodes_j) = Kg(enodes_i, enodes_j) + kConv(ii, j);
				K_gFlux_inlet.coeffRef(enodes_i, enodes_j) +=  + kConv(ii, j); //Uncomment 
    		}
    		
		}

		dK_gFlux_inlet_drho1[el] += dkConv_drho1;
		dK_gFlux_inlet_drho2[el] += dkConv_drho2;
	}

	VectorXi outletEl = x_y_front;
	for (int i = 0; i < outletEl.rows(); i++) {
		int el = outletEl(i);
		MatrixXd enodes = nodemap.row(el);
		MatrixXd elem_coords(enodes.cols(), dim); // convert dim to int, dim=3
		for (int ii = 0; ii < enodes.cols(); ii++) {
        	elem_coords.row(ii) = coordinates.row(enodes(ii));
			int enodes_i = enodes(ii);
			L_outlet(enodes_i) = 1;
    	}

		double w = (1 - Rho(el,1))*Rho(el,0);
		double dw_drho1 = 1 - Rho(el,1);
		double dw_drho2 = - Rho(el,0);

		MatrixXd wtEdge(1, 4);
		wtEdge << 1, 1, 1, 1; //Gauss Quadrature Weights

		MatrixXd kConv;
		kConv = MatrixXd::Zero(pow(2,dim), pow(2,dim)); //Change 8 to 2**dim
		MatrixXd dkConv_drho1;
		dkConv_drho1 = MatrixXd::Zero(pow(2,dim), pow(2,dim));
		MatrixXd dkConv_drho2;
		dkConv_drho2 = MatrixXd::Zero(pow(2,dim), pow(2,dim));

		MatrixXd gPointsEdge(4, dim);

		gPointsEdge << -1/sqrt(3), -1/sqrt(3), -1,
               			1/sqrt(3), -1/sqrt(3), -1,
               			1/sqrt(3), 1/sqrt(3), -1,
               	   		-1/sqrt(3), 1/sqrt(3), -1;

		for (int aa = 0; aa < pow(2,dim-1); aa++) {
			MatrixXd N;
			MatrixXd dN;
			shapeFcn(aa, gPointsEdge, myMesh, N, dN);
			//std::cout << "ipt_coords: \n" << ipt_coords << std::endl;
			MatrixXd J = elem_coords.transpose() * dN;
			double dS = sqrt(pow(J(0,0),2) + pow(J(1,0),2)); //Check dS for 3D
			kConv = kConv + w * N * N.transpose() * dS * wtEdge(aa); //where hConv = 1
			dkConv_drho1 = dkConv_drho1 + dw_drho1 * N * N.transpose() * dS * wtEdge(aa);
			dkConv_drho2 = dkConv_drho2 + dw_drho2 * N * N.transpose() * dS * wtEdge(aa);

		}

		for (int ii = 0; ii < enodes.cols(); ++ii) {
			int enodes_i = enodes(ii);
    		for (int j = 0; j < enodes.cols(); ++j) {
				int enodes_j = enodes(j);
				//Kg(enodes_i, enodes_j) = Kg(enodes_i, enodes_j) + kConv(ii, j);
				K_gFlux_outlet.coeffRef(enodes_i, enodes_j) +=  + kConv(ii, j); //Uncomment 
    		}
    		
		}

		dK_gFlux_outlet_drho1[el] += dkConv_drho1;
		dK_gFlux_outlet_drho2[el] += dkConv_drho2;
	}

	MatrixXd f2 = (L_inlet.transpose() * K_gFlux_inlet * T) - (L_outlet.transpose() * K_gFlux_outlet * T);
	std::cout << "HT_inlet: " << (L_inlet.transpose() * K_gFlux_inlet * T) << std::endl;
	std::cout << "HT_outlet: " << (L_outlet.transpose() * K_gFlux_outlet * T) << std::endl;

	MatrixXd dfdT = (L_inlet.transpose() * K_gFlux_inlet) - (L_outlet.transpose() * K_gFlux_outlet);
	Eigen::SparseMatrix<double> dR_thf_dTf = myGrad_Helper_Vars.Kth_freefree;
	MatrixXd dfdTf = Eigen::MatrixXd::Zero(th_freedofs.cols(), 1);

	for (int c = 0; c < th_freedofs.cols(); c++) {
		dfdTf(c) = dfdT(th_freedofs(c));
	}

	
	MatrixXd AdjTh = Eigen::MatrixXd::Zero(ndof, 1);
	BiCGSTAB<SparseMatrix<double> > solver;
	solver.compute(dR_thf_dTf);
	MatrixXd AdjTh_f = solver.solve(-dfdTf);
	
	for (int c = 0; c < th_freedofs.cols(); c++) {
		AdjTh(th_freedofs(c)) = AdjTh_f(c);
	}

	myGrad_Helper_Vars.df2drho1.resize(1, nel);
	myGrad_Helper_Vars.df2drho2.resize(1, nel);

	for (int i = 0; i < nel; ++i) {
		MatrixXd enodes = nodemap.row(i);
		
		MatrixXd T_el(enodes.cols(), 1); 
		MatrixXd Adjth_el(enodes.cols(), 1); 
		MatrixXd L_inlet_el(enodes.cols(), 1);
		MatrixXd L_outlet_el(enodes.cols(), 1);
		for (int ii = 0; ii < enodes.cols(); ii++) {
			int enodes_i = enodes(ii);
        	T_el(ii) = T(enodes_i);
			Adjth_el(ii) = AdjTh(enodes_i);
			L_inlet_el(ii) = L_inlet(enodes_i);
			L_outlet_el(ii) = L_outlet(enodes_i);
    	}
		
		MatrixXd delfdelrho1 = (L_inlet_el.transpose() * dK_gFlux_inlet_drho1[i] * T_el) - (L_outlet_el.transpose() * dK_gFlux_outlet_drho1[i] * T_el);
		MatrixXd delfdelrho2 = (L_inlet_el.transpose() * dK_gFlux_inlet_drho2[i] * T_el) - (L_outlet_el.transpose() * dK_gFlux_outlet_drho2[i] * T_el);

		myGrad_Helper_Vars.df2drho1(i) = delfdelrho1(0,0) + (Adjth_el.transpose()*(myGrad_Helper_Vars.dKthdrho1[i]*T_el - myGrad_Helper_Vars.dPthdrho1[i]))(0,0);
		myGrad_Helper_Vars.df2drho2(i) = delfdelrho2(0,0) + (Adjth_el.transpose()*(myGrad_Helper_Vars.dKthdrho2[i]*T_el - myGrad_Helper_Vars.dPthdrho2[i]))(0,0);

		
	}

	return f2(0,0);
	// return 0;

}

std::vector<std::pair<std::string, std::vector<int>>> read_csv2(std::string filename){
    // Reads a CSV file into a vector of <string, vector<int>> pairs where
    // each pair represents <column name, column values>

    // Create a vector of <string, int vector> pairs to store the result
    std::vector<std::pair<std::string, std::vector<int>>> result;

    // Create an input filestream
    std::ifstream myFile(filename);

    // Make sure the file is open
    if(!myFile.is_open()) throw std::runtime_error("Could not open file");

    // Helper vars
    std::string line, colname;
    int val;

    // Read the column names
    /*if(myFile.good())
    {
        // Extract the first line in the file
        std::getline(myFile, line);

        // Create a stringstream from line
        std::stringstream ss(line);

        // Extract each column name
        while(std::getline(ss, colname, ',')){
            
            // Initialize and add <colname, int vector> pairs to result
            result.push_back({colname, std::vector<int> {}});
        }
    }*/

    // Read data, line by line
    while(std::getline(myFile, line))
    {
        // Create a stringstream of the current line
        std::stringstream ss(line);
        
        // Keep track of the current column index
        int colIdx = 0;
        
        // Extract each integer
        while(ss >> val){
            
            // Add the current integer to the 'colIdx' column's values vector
            result.at(colIdx).second.push_back(val);
            
            // If the next token is a comma, ignore it and move on
            if(ss.peek() == ',') ss.ignore();
            
            // Increment the column index
            colIdx++;
        }
    }

    // Close file
    myFile.close();

    return result;
}

std::vector<std::vector<double>> read_csv(const std::string& filename) {
    // Create a vector of vectors to store the result
    std::vector<std::vector<double>> result;

    // Create an input filestream
    std::ifstream myFile(filename);

    // Make sure the file is open
    if (!myFile.is_open()) {
        std::cerr << "Could not open file: " << filename << std::endl;
        return result;
    }

    // Helper vars
    std::string line;
    double val;

    // Read data, line by line
    while (std::getline(myFile, line))
    {
        // Create a stringstream of the current line
        std::stringstream ss(line);

        // Create a vector to store the values of the current line
        std::vector<double> row;

        // Extract each double
        while (ss >> val) {
            row.push_back(val);

            // If the next token is a comma, ignore it and move on
            if (ss.peek() == ',') ss.ignore();
        }

        // Add the row to the result vector
        result.push_back(row);
    }

    // Close file
    myFile.close();

    return result;
}