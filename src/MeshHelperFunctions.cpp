#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include "MeshHelperFunctions.h"
#include <chrono>

using namespace Eigen;


Eigen::MatrixXd generateGridCoordinates(double L, double H, double B, int nelx, int nely, int nelz) {
    // Create Eigen vectors for x, y, and z
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(nelx + 1, 0, L);
    Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(nely + 1, 0, H);
    Eigen::VectorXd z = Eigen::VectorXd::LinSpaced(nelz + 1, 0, B);

	MatrixXd coordinates;
	coordinates = MatrixXd::Zero((nelx + 1) * (nely + 1) * (nelz + 1), 3);

    // Output the coordinates
	int iter = 0;

    for (int k = 0; k <= nelz; ++k) {
		//for (int j = 0; j <= nely; ++j) {
        for (int i = 0; i <= nelx; ++i) {
            for (int j = nely; j >= 0; --j) {
                //std::cout << "[" << x[i] << ", " << y[j] << ", " << z[k] << "],\n";
				coordinates(iter,0) = x(i);
				coordinates(iter,1) = y(j);
				coordinates(iter,2) = z(k);
				 iter ++;
            }
        }
    }

	return coordinates;
}

Eigen::MatrixXd generateNodeMap(int nelx, int nely, int nelz) {
    // Create Node Map Matrix
	MatrixXd nodemap;
	nodemap = MatrixXd::Zero(nelx * nely * nelz, 8);
	int nnodex = nelx + 1;
	int nnodey = nely + 1;
	int nnodez = nelz + 1;

    // Assign Node Map
	int iter = 0;

	for (int k = 0; k < nelz; k++) {
        for (int i = 0; i < nelx; i++) {
            for (int j = 0; j < nely; j++) {
                //std::cout << "[" << x[i] << ", " << y[j] << ", " << z[k] << "],\n";
				nodemap(iter,0) = nnodey*i + (j+1) + (nnodex*nnodey*k);
				nodemap(iter,1) = nnodey*(i+1) + (j+1) + (nnodex*nnodey*k);
				nodemap(iter,2) = nnodey*(i+1) + (j) + (nnodex*nnodey*k);
				nodemap(iter,3) = nnodey*i + (j) + (nnodex*nnodey*k);
				nodemap(iter,4) = nnodey*i + (j+1) + (nnodex*nnodey*(k+1));
				nodemap(iter,5) = nnodey*(i+1) + (j+1) + (nnodex*nnodey*(k+1));
				nodemap(iter,6) = nnodey*(i+1) + (j) + (nnodex*nnodey*(k+1));
				nodemap(iter,7) = nnodey*i + (j) + (nnodex*nnodey*(k+1));
				 iter ++;
            }
        }
    }


	return nodemap;
}

Eigen::MatrixXd generateMeshCenters(const Eigen::MatrixXd nodemap, const Eigen::MatrixXd coordinates, int nel) {
    // Create Node Map Matrix
	MatrixXd mesh_centers;
	mesh_centers = MatrixXd::Zero(nel, 3);

	for (int i = 0; i < nel; i++) {
		int node0 = nodemap(i, 0);
		int node1 = nodemap(i, 1);
		int node3 = nodemap(i, 3);
		int node4 = nodemap(i, 4);
		mesh_centers(i,0) = 0.5 * (coordinates(node0,0) + coordinates(node1,0));
		mesh_centers(i,1) = 0.5 * (coordinates(node0,1) + coordinates(node3,1));
		mesh_centers(i,2) = 0.5 * (coordinates(node0,2) + coordinates(node4,2));
		if (i < 51) {
			/*std::cout << "coordinates(node0,0): " << coordinates(node0,0) << std::endl;
			std::cout << "coordinates(node1,0): " << coordinates(node1,0) << std::endl;
			std::cout << "coordinates(node0,1): " << coordinates(node0,1) << std::endl;
			std::cout << "coordinates(node3,1): " << coordinates(node3,1) << std::endl;
			std::cout << "coordinates(node0,2): " << coordinates(node0,2) << std::endl;
			std::cout << "coordinates(node4,2): " << coordinates(node4,2) << std::endl;*/
			std::cout << "Element " << i << ": (" << mesh_centers(i,0) << " " << mesh_centers(i,1) << " " << mesh_centers(i,2) << ")" << std::endl;
		}
    }


	return mesh_centers;
}

Eigen::SparseMatrix<double> generateElemToNodalDensityMatrix(const Eigen::MatrixXd nodemap, int nel) {
// Eigen::MatrixXd generateElemToNodalDensityMatrix(const Eigen::MatrixXd nodemap, int nel) {
    // Create Node Map Matrix
	std::cout << "generateElemToNodalDensityMatrix" <<std::endl;
	auto start_time = std::chrono::high_resolution_clock::now();
	int nnodes = nodemap.maxCoeff() + 1;
	// MatrixXd H; H = MatrixXd::Zero(nnodes, nel); // Dense
	Eigen::SparseMatrix<double> H(nnodes, nel); //Sparse

	for (int i = 0; i < nel; i++) {
		MatrixXd enodes = nodemap.row(i);
		//std::cout << "enodes: " << enodes << std::endl;
		for (int j = 0; j < enodes.cols(); j++) {
			int enodes_j = enodes(j);
			H.insert(enodes_j,i) = 1; //Sparse
			// H.(enodes_j,i) = 1; //Dense 
    	}
    }
	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "Creating H time: " << duration.count() << " seconds" << std::endl;


	// Dense Matrix Implementation
	/*VectorXd H_rowSums = H.rowwise().sum();
	std::cout << "H:\n" << H <<std::endl;
	std::cout << "H_rowSums:\n" << H_rowSums <<std::endl;
	MatrixXd Havg = H.array().colwise() / H_rowSums.array();*/

	
	// Sparse Implementation
	start_time = std::chrono::high_resolution_clock::now();
	Eigen::SparseMatrix<double> Havg(nnodes, nel);
	Havg.reserve(H.nonZeros()); // Reserve space for non-zero elements
	std::cout << "H.nonZeros(): " << H.nonZeros() << std::endl;

	for (int i = 0; i < nnodes; ++i) {
        // Extract the column and find the non-zero indices
        Eigen::SparseVector<double> row = H.row(i);
		//std::cout << "H.row(" << i << "):" << row << " row.sum(): " << row.sum() << std::endl;
        for (Eigen::SparseVector<double>::InnerIterator it(row); it; ++it) {
            // Compute the element-wise quotient and assign it to Havg
            Havg.insert(i, it.row()) = it.value() / row.sum();
			//std::cout << "it: " << it.value() << " it.col(): " << it.col() << " Division: " << it.value() / row.sum() << std::endl;
        }
    }
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "Creating Havg time: " << duration.count() << " seconds" << std::endl;


	// std::cout << "Havg(19,5):\n" << Havg.coeff(19,5) <<std::endl;

	return Havg;

	// Alternate Slower Sparse Matrix Implementation
	/*Eigen::SparseVector<double> H_rowSums(nnodes);
	for (int i = 0; i < nel; ++i) {
        H_rowSums += H.col(i);
    }
	std::cout << "H_rowSums:\n" << H_rowSums <<std::endl;
	// Check if you can do this without Loops
	Eigen::SparseMatrix<double> Havg(nnodes, nel);
	Havg.reserve(H.nonZeros()); // Reserve space for non-zero elements
	for (int i = 0; i < nel; ++i) {
		Havg.col(i) = H.col(i).cwiseQuotient(H_rowSums);
	}*/
}