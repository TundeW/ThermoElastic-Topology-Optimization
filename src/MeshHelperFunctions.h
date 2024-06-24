#pragma once
#ifndef MESHHELPERFUNCTIONS_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define MESHHELPERFUNCTIONS_H

Eigen::MatrixXd generateGridCoordinates(double L, double H, double B, int nelx, int nely, int nelz);
Eigen::MatrixXd generateNodeMap(int nelx, int nely, int nelz);
Eigen::MatrixXd generateMeshCenters(const Eigen::MatrixXd nodemap, const Eigen::MatrixXd coordinates, int nel);
// Eigen::MatrixXd generateElemToNodalDensityMatrix(const Eigen::MatrixXd nodemap, int nel);
Eigen::SparseMatrix<double> generateElemToNodalDensityMatrix(const Eigen::MatrixXd nodemap, int nel);

#endif