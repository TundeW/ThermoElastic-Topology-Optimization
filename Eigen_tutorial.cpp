/*
#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

int main() {
	std::cout << "Hello World!" << std::endl;
	std::cin.get();
	Matrix <float, 3,3> matrixA;
	matrixA.setZero();

	std::cout << matrixA << std::endl;

	std::cout << "Changing matrix entries" << std::endl;

	matrixA(0,0) = 2;
	matrixA(0,1) = 3;
	matrixA(2,1) = 5;

	std::cout << matrixA << std::endl;

	// Define a 3x3 matrix of floats and set its entries to zero -typedef declaration
	//Used for small sized matrices, whose size we know before compile time
	Matrix3f matrixA1;
	matrixA1.setZero();
	std::cout << matrixA1 << std::endl;

	// Define a dynamic matrix, explicit declaration
	Matrix <float, Dynamic, Dynamic> matrixB;

	//Define a dynamic matrix, typedef declaration
	MatrixXf matrixB1;

	//Constructor for allocating memory but not initializing it
	MatrixXd matrixC(10,10);

	// Accessing and printing entries of a matrix
	Matrix4f matrixD;
	matrixD << 1,   2,  3,  4,
    	    5,   6,  7,  8,
        	9,  10, 11, 12,
        	13, 14, 15, 16;
	std::cout << std::endl <<  matrixD << std::endl;

	// define a dynamic matrix, resize it to a 3x3 matrix, 
	// and set its entries to zero, and print the matrix
	Matrix3f matrixD1;
	matrixD1.resize(3, 3);
	matrixD1.setZero(3, 3);
	std::cout << std::endl << matrixD1;

	// setting matrix entries - two approaches
    int rowsNumber = 5;
    int columnNumber= 5;

	// matrix of zeros
    MatrixXf matrix1zeros;
	matrix1zeros = MatrixXf::Zero(rowsNumber, columnNumber);
    std::cout << "\n \n" << matrix1zeros << std::endl;

	// another option:
    MatrixXf matrix1zeros1;
    matrix1zeros1.setZero(rowsNumber, columnNumber);
    std::cout << "\n \n" << matrix1zeros1 << std::endl;

	//matrix of ones
    MatrixXf matrix1ones;
    matrix1ones = MatrixXf::Ones(rowsNumber, columnNumber);
    std::cout << "\n \n" << matrix1ones << std::endl;

    //another option
    MatrixXf matrix1ones1;
    matrix1ones1.setOnes(rowsNumber, columnNumber);
    std::cout << "\n \n" << matrix1ones1 << std::endl;

	//matrix of constants
    float value = 1.1;
    MatrixXf matrix1const;
    matrix1const = MatrixXf::Constant(rowsNumber, columnNumber,value);
    std::cout << "\n \n" << matrix1const << std::endl;

    //another option
    MatrixXf matrix1const1;
    matrix1const1.setConstant(rowsNumber, columnNumber, value);
    std::cout << "\n \n" << matrix1const1 << std::endl;

	//identity matrix, two approaches
    int rowNumber = 10;
    columnNumber = 10;
     
    MatrixXd matrixIdentity;
    matrixIdentity = MatrixXd::Identity(rowNumber,columnNumber);
    std::cout << "\n \n" << matrixIdentity << std::endl;

	MatrixXd matrixIdentity1;
    matrixIdentity1.setIdentity(rowNumber, columnNumber);
    std::cout << "\n \n" << matrixIdentity1 << std::endl;

	//accessing matrix blocks
    MatrixXd matrixV(4,4);
    matrixV << 101, 102, 103, 104,
    	 	   105, 106, 107, 108,
    		   109, 110, 111, 112,
     		   113, 114, 115, 116;
    //access the matrix composed of 1:2 rows and 1:2 columns of matrixV
    MatrixXd matrixVpartition = matrixV.block(0, 0, 2, 2);
    std::cout << "\n \n" << matrixVpartition << std::endl;
     
    MatrixXd matrixVpartition2 = matrixV.block(1,1, 2, 2);
    std::cout << "\n \n" << matrixVpartition2 << std::endl;

	//accessing columns and rows of a matrix
 
	std::cout << "\n\n"<<"Row 1 of matrixV is \n " << matrixV.row(0);
	std::cout << "\n\n" << "Column 1 of matrixV is \n" << matrixV.col(0);
	
	//create a diagonal matrix out of a vector
    Matrix <double, 3, 1> vector1;
    vector1 << 1, 2, 3;
    MatrixXd diagonalMatrix;
    diagonalMatrix = vector1.asDiagonal();
    std::cout << " Diagonal matrix is \n\n" << diagonalMatrix;

	// basic matrix operations
    //matrix addition
    MatrixXd A1(2, 2);
    MatrixXd B1(2, 2);
 
    A1 << 1, 2,
        3, 4;
    B1 << 3, 4,
        5, 6;
    MatrixXd C1 = A1 + B1;
    std::cout << " \n\n The sum of A1 and B1 is\n\n" << C1 << std::endl;
    // similarly you can subtract A1 and B1

	// Matrix Multiplication
    MatrixXd D1 = A1 * B1;
    std::cout << " \n\n The matrix product of A1 and B1 is\n\n" << D1 << std::endl;
 
    //multiplication by a scalar
    int s1 = 2;
    MatrixXd F1;
    F1 = s1 * A1;
    std::cout << " \n\n The matrix product of the scalar 2 and the matrix A1 is\n\n" << F1 << std::endl;
	
	//Matrix Transpose
    MatrixXd At;
    MatrixXd R1;
    // we can obain a transpose of A1 as follows
    At = A1.transpose();
    std::cout << "\n\n Original matrix A1\n\n" << A1;
    std::cout << "\n\n Its transpose\n\n " << At;
 
    // we can use a transpose operator in expressions
    R1 = A1.transpose() + B1;
    std::cout << "\n\n Matrix R1=A1^{T}+B1 is\n\n" << R1;
 
    // here we should be careful, the operation .transpose() might lead to unexpected results in this scenarios
    //  A1 = A1.transpose();
    // std::cout << " \n\n This should be a transpose of the matrix A1\n\n" << A1 << std::endl;
 
    // the correct and safe way to do the matrix transpose is the following
    A1.transposeInPlace();
    std::cout << " \n\n This should be a transpose of the matrix A1\n\n" << A1 << std::endl;
    //restore A1 matrix to its original state
    A1.transposeInPlace();

	// Matrix Inverse
	MatrixXd G1;
	G1 = A1.inverse();
	std::cout << "\n\n The inverse of the matrix A1 is\n\n" << G1;
	std::cout << "\n\n Double check A1^{-1}*A1=\n\n" << G1*A1;
}
*/

/*void generateGridCoordinates(double L, double H, double B, int nelx, int nely, int nelz) {
    // Create Eigen vectors for x, y, and z
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(nelx + 1, 0, L);
    Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(nely + 1, 0, H);
    Eigen::VectorXd z = Eigen::VectorXd::LinSpaced(nelz + 1, 0, B);

    // Create meshgrid using Eigen's broadcasting
    Eigen::MatrixXd X = x.replicate(1, nely + 1);
    Eigen::MatrixXd Y = y.replicate(nelx + 1, 1);
    Eigen::MatrixXd Z = z.replicate(nelx + 1, nely + 1).transpose();

    // Combine X, Y, Z to get the final coordinates
    Eigen::MatrixXd coordinates(X.rows() * Y.cols(), 3);
    coordinates << X.array().flatten(), Y.array().flatten(), Z.array().flatten();

    // Output the coordinates
    std::cout << "Coordinates:\n" << coordinates << "\n";
}*/