#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include "ProjectionHelperFunctions.h"

using namespace Eigen;

/*typedef struct {
    Eigen::MatrixXd coordinate;
	Eigen::MatrixXd nodemap;
	Eigen::MatrixXd mesh_centers;
	double nel;
	int dim;
	double L;
	double H;
	double B;
	double nelx;
	double nely;
	double nelz;
	double Area;
} MESH;*/

Eigen::MatrixXd BarProjection(const Eigen::MatrixXd Centre, const Eigen::VectorXd R, const MESH myMesh, Eigen::MatrixXd& dRhoodx, Eigen::MatrixXd& dRhoodr, Eigen::MatrixXd& dRhood_) {
    // Access struct parameters
    Eigen::MatrixXd coordinates = myMesh.coordinate;
    Eigen::MatrixXd nodemap = myMesh.nodemap;
    Eigen::MatrixXd mesh_centers = myMesh.mesh_centers;
	double nel = myMesh.nel;

	//std::cout << "Coordinates:\n";
	//std::cout << coordinates << std::endl;
	//std::cout << "Nodemap:\n" << nodemap << "\n";
	//std::cout << "Mesh Centers:\n" << mesh_centers << "\n";
	//std::cout << "nel: " << nel << std::endl;

	// Create Node Map Matrix
	MatrixXd Rho;
	Rho = MatrixXd::Zero(nel, 1);

	double rad = 0.02;
	rad = 0.003125;
	rad = 0.0025*1.25;
	//rad = 0.15;
	double rho_min = 1e-3; 

	MatrixXd D;
	D = MatrixXd::Zero(mesh_centers.rows(), Centre.rows());

	MatrixXd Phi_q;
	Phi_q = MatrixXd::Zero(mesh_centers.rows(), Centre.rows());
	MatrixXd DPhi_qDx;
	DPhi_qDx = MatrixXd::Zero(mesh_centers.rows(), 2*Centre.rows());

	MatrixXd DPhi_qDr;
	DPhi_qDr = MatrixXd::Zero(mesh_centers.rows(), Centre.rows());

	for (int i = 0; i < Centre.rows(); i++) {
		MatrixXd centre = Centre.row(i);
		double r = R(i);
		//std::cout << "Pipe " << i << " Centre: " << centre << " Radius: " << r << std::endl;
		//MatrixXd x0 = centre;
		//MatrixXd xv = centre;

		 // Append 0 to each row
    	MatrixXd x0(centre.rows(), centre.cols() + 1);
    	//x0 << centre, Eigen::MatrixXd::Zero(centre.rows(), 1);
		x0 << centre, 0;

    	// Append 2 to each row
    	MatrixXd xf(centre.rows(), centre.cols() + 1);
    	//xf << centre, Eigen::MatrixXd::Constant(centre.rows(), 1, 2.0);
		xf << centre, 2;

		// Repeat x0 and xffor each row in mesh_centers
    	MatrixXd x0v = x0.replicate(mesh_centers.rows(), 1);
		MatrixXd xfv = xf.replicate(mesh_centers.rows(), 1);

		// Preprocess
		MatrixXd a_b = xfv - x0v;
		MatrixXd l_b = a_b.rowwise().norm();
		MatrixXd a_bar = a_b.rowwise().normalized();

		MatrixXd b = mesh_centers - x0v;
		MatrixXd e = mesh_centers - xfv;
		MatrixXd norm_x_e_1b = b.rowwise().norm();
		MatrixXd norm_x_e_2b = e.rowwise().norm();

		MatrixXd b_dot_a_bar = b.array() * a_bar.array();
		MatrixXd l_be = b_dot_a_bar.rowwise().sum();

		MatrixXd l_be_dot_a_bar(l_be.rows(), a_bar.cols());
		for (int i = 0; i < l_be.rows(); ++i) {
    		l_be_dot_a_bar.row(i) = l_be(i) * a_bar.row(i).array();
		}

    	// Print the result
    	//std::cout << "l_be_dot_a_bar:\n" << l_be_dot_a_bar << std::endl;
		//std::cout << "l_b:\n" << l_b << std::endl;
		//std::cout << "l_be:\n" << l_be << std::endl;
		MatrixXd vec_r_be  = b - l_be_dot_a_bar;
		MatrixXd r_be = vec_r_be.rowwise().norm();

		//ArrayXd branch1 = l_be.array() <= 0.0;
		//MatrixXd branch1 = l_be.cast<double>().array() <= 0.0;
		//MatrixXd branch1 = (l_be.cast<double>().array() <= 0.0).matrix();
		//ArrayXXf branch1 = (l_be.array() <= 0.0).cast<float>()*2.f-1.f;
		//ArrayXXf branch2 = (l_be.array() > l_b.array()).cast<float>()*2.f-1.f;
		//ArrayXXf branch3 = (!(branch1 || branch2)).cast<float>()*2.f-1.f;
		//ArrayXd branch2 = l_be.array() <= l_b.array();
		//ArrayXd branch3 = !(branch1 || branch2);

		MatrixXd dd;
		dd = MatrixXd::Zero(mesh_centers.rows(), 1);

		MatrixXd phi_q;
		phi_q = MatrixXd::Zero(mesh_centers.rows(), 1);

		MatrixXd dPhi_qdx0;
		dPhi_qdx0 = MatrixXd::Zero(mesh_centers.rows(), 3);

		MatrixXd dPhi_qdxf;
		dPhi_qdxf = MatrixXd::Zero(mesh_centers.rows(), 3);

		MatrixXd dPhi_qdx;
		dPhi_qdx = MatrixXd::Zero(mesh_centers.rows(), 2);

		MatrixXd dPhi_qdr;
		dPhi_qdr = -1 * MatrixXd::Constant(mesh_centers.rows(), 1, 1.0);


		//std::cout << "b:\n" << b <<std::endl; 
		//std::cout << "b2:\n" << -b.row(0)/norm_x_e_1b(0) <<std::endl; 

		for (int j = 0; j < l_be.rows(); j++) {
			if (l_be(j) <= 0) {
				dd(j) = norm_x_e_1b(j);
				dPhi_qdx0.row(j) = -b.row(j)/norm_x_e_1b(j);
				//std::cout << "norm_x_e_1b" << std::endl;
			}
			else if (l_be(j) > l_b(j)) {
				dd(j) = norm_x_e_2b(j);
				dPhi_qdxf.row(j) = -e.row(j)/norm_x_e_2b(j);
				//std::cout << "norm_x_e_2b" << std::endl;
			}
			else {
				dd(j) = r_be(j);
				/*dPhi_qdx0(j,0) = - (1 + l_be(j)/l_b(j)*(-1));
				dPhi_qdx0(j,1) = - (1 + l_be(j)/l_b(j)*(-1));
				dPhi_qdx0(j,2) = - (1 + l_be(j)/l_b(j)*(-1));
				dPhi_qdxf(j,0) = - (0 + l_be(j)/l_b(j)*(1));
				dPhi_qdxf(j,1) = - (0 + l_be(j)/l_b(j)*(1));
				dPhi_qdxf(j,2) = - (0 + l_be(j)/l_b(j)*(1));*/

				dPhi_qdx0.row(j) = (-vec_r_be.row(j)/r_be(j)) * (1 + l_be(j)/l_b(j)*(-1));
				dPhi_qdxf.row(j) = (-vec_r_be.row(j)/r_be(j)) * (0 + l_be(j)/l_b(j)*(1));
				
				//std::cout << "r_be" << std::endl;
			}

			phi_q(j) = dd(j) - r;
		}
		
		MatrixXd dx0_fdx;
		dx0_fdx.resize(3, 2);
    	dx0_fdx << 1, 0,
              	   0, 1,
              	   0, 0;
		
		dPhi_qdx = dPhi_qdx0*dx0_fdx + dPhi_qdxf*dx0_fdx;
		//std::cout << "dPhi_qdx0:\n" << dPhi_qdx0.transpose() << std::endl;
		//std::cout << "dPhi_qdxf:\n" << dPhi_qdxf.transpose() << std::endl;

		D.col(i) = dd;
		Phi_q.col(i) = phi_q;
		DPhi_qDx.col(2 * i) = dPhi_qdx.col(0);
		DPhi_qDx.col(2 * i + 1) = dPhi_qdx.col(1);
		DPhi_qDr.col(i) = dPhi_qdr;
	}

	//std::cout << "DPhi_qDx:\n" << DPhi_qDx << std::endl;

	MatrixXd dPhi_qdx;
	dPhi_qdx = MatrixXd::Zero(mesh_centers.rows(), 2*Centre.rows());

	MatrixXd dPhi_qdr;
	dPhi_qdr = MatrixXd::Zero(mesh_centers.rows(), Centre.rows());

	VectorXd d = D.rowwise().minCoeff();
	//VectorXd phi_q = Phi_q.rowwise().minCoeff();
	VectorXd phi_q(Phi_q.rows());
	VectorXi pipenum(phi_q.rows());
  	for (int i = 0; i < phi_q.rows(); ++i) {
		phi_q(i) = Phi_q.row(i).minCoeff( &pipenum[i] );

		dPhi_qdx(i,2*pipenum(i)) = DPhi_qDx(i,2*pipenum(i));
		dPhi_qdx(i,2*pipenum(i)+1) = DPhi_qDx(i,2*pipenum(i)+1);
		dPhi_qdr(i,pipenum(i)) = DPhi_qDr(i,pipenum(i));

  	}
	//std::cout << "pipenum:\n" << pipenum << std::endl;

	

	//std::cout << "Phi_q:\n" << Phi_q << std::endl;
	//std::cout << "phi_q:\n" << phi_q << std::endl;

	//std::cout << "dPhi_qdx:\n" << dPhi_qdx << std::endl;

	/*Eigen::VectorXi pipenum = Eigen::VectorXi::Zero(Phi_q.rows());

    for (int i = 0; i < Phi_q.rows(); ++i) {
        double minVal;
        int minIndex;
        Phi_q.row(i).minCoeff(&minVal, &minIndex);
        pipenum(i) = minIndex;
    }*/
	MatrixXd rho_q;
	rho_q = MatrixXd::Zero(mesh_centers.rows(), 1);

	MatrixXd drho_dd;
	drho_dd = MatrixXd::Zero(mesh_centers.rows(), 1);

	MatrixXd dRhodx;
	dRhodx = MatrixXd::Zero(mesh_centers.rows(), 2*Centre.rows());

	MatrixXd dRhodr;
	dRhodr = MatrixXd::Zero(mesh_centers.rows(), Centre.rows());

	//MatrixXd Rho;
	//Rho = MatrixXd::Zero(mesh_centers.rows(), 1);

	for (int i = 0; i < mesh_centers.rows(); i++) {
		if (phi_q(i) > rad) {
			rho_q(i) = 0;
		}
		else if (abs(phi_q(i)) <= rad) {
			//rho_q(i) = (1/(np.pi*rad**2))*(rad**2*np.arccos(phi_q[idxR2]/rad)-phi_q[idxR2]*np.sqrt(rad**2-phi_q[idxR2]**2));
			rho_q(i) = (1 / (M_PI * rad * rad)) * (rad * rad * acos(phi_q(i) / rad) - phi_q(i) * sqrt(rad * rad - phi_q(i) * phi_q(i)));
			drho_dd(i) = (-2/(M_PI * rad * rad)) * sqrt(rad * rad - phi_q(i) * phi_q(i));
		}
		else if (phi_q(i) < -rad) {
			rho_q(i) = 1;
		}

		dRhodx.row(i) = drho_dd(i)*dPhi_qdx.row(i);
		dRhodr.row(i) = drho_dd(i)*dPhi_qdr.row(i);
		Rho(i) = rho_min + rho_q(i)*(1-rho_min);
	}
	// Print the results
    //std::cout << "Original Centre:\n" << centre << " " << centre.rows() << " " << centre.cols() << "\n\n";
    //std::cout << "phi_q:\n" << phi_q << "\n";
	//std::cout << " branch1:\n" << branch1 << "\n";
	//std::cout << " branch2:\n" << branch2 << "\n";
	//std::cout << " branch3:\n" << branch3 << "\n";
	//std::cout << " dRhodx:\n" << dRhodx << "\n";
	//std::cout << " dRhodr:\n" << dRhodr << "\n";

	//MatrixXd dRhod_;
	dRhoodx = dRhodx;
	dRhoodr = dRhodr;
	dRhood_ = MatrixXd::Zero(mesh_centers.rows(), Centre.rows());

	return rho_q;
}

Eigen::MatrixXd PipeProjection(const Eigen::MatrixXd Centre, const Eigen::VectorXd Inner_r, const Eigen::VectorXd Outer_r, const MESH myMesh, GRAD_HELPER_VARS& myGrad_Helper_Vars) {
    // Access struct parameters
    Eigen::MatrixXd coordinates = myMesh.coordinate;
    Eigen::MatrixXd nodemap = myMesh.nodemap;
    Eigen::MatrixXd mesh_centers = myMesh.mesh_centers;
	double nel = myMesh.nel;

	// Create Node Map Matrix
	MatrixXd Rho;
	Rho = MatrixXd::Zero(nel, 2);

	MatrixXd dinner_rho_dc;
	MatrixXd dinner_rho_dr;
	MatrixXd dinner_rho_dt;
	MatrixXd douter_rho_dc;
	MatrixXd douter_rho_dr;
	MatrixXd douter_rho_dt;

	double rho_min = 0;
	MatrixXd inner_rho = BarProjection(Centre, Inner_r, myMesh, dinner_rho_dc, dinner_rho_dr, dinner_rho_dt);
	MatrixXd outer_rho = BarProjection(Centre, Outer_r, myMesh, douter_rho_dc, douter_rho_dt, douter_rho_dr);

	Eigen::MatrixXd dinner_rho_dx(dinner_rho_dc.rows(), dinner_rho_dc.cols() + dinner_rho_dr.cols() + dinner_rho_dt.cols());
	dinner_rho_dx << dinner_rho_dc, dinner_rho_dr, dinner_rho_dt;

	Eigen::MatrixXd douter_rho_dx(douter_rho_dc.rows(), douter_rho_dc.cols() + douter_rho_dr.cols() + douter_rho_dt.cols());
	douter_rho_dx << douter_rho_dc, douter_rho_dr, douter_rho_dt;

	//std::cout << "dinner_rho_dx:\n" << dinner_rho_dx << std::endl;
	//std::cout << "douter_rho_dx:\n" << douter_rho_dx << std::endl;

	MatrixXd rho1 = 1 - (outer_rho - inner_rho).array();
	MatrixXd rho2 = 1 - inner_rho.array();
	//MatrixXd rho2 = outer_rho.array();
	//rho2 = 1 - outer_rho.array();

	// Convert to density
	MatrixXd Rho1 = rho_min + (1-rho_min)*rho1.array();
    MatrixXd Rho2 = rho_min + (1-rho_min)*rho2.array();

	myGrad_Helper_Vars.dRho1_dx = (1-rho_min)*(-(douter_rho_dx - dinner_rho_dx));
	myGrad_Helper_Vars.dRho2_dx = (1-rho_min)*(-dinner_rho_dx);

	//std::cout << "Rho1:\n" << Rho1 << "\n\n";
	//std::cout << "Rho2:\n" << Rho2 << "\n\n";
	
	//Rho.col(0) = Rho1;
	//Rho.col(1) = Rho2;
	Rho.col(0) = Rho1;
	Rho.col(1) = Rho2;

	return Rho;
}