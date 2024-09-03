#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <fstream>
#include "WriteOutputFunctions.h"

using namespace Eigen;

void WriteObjective(std::ofstream &objectiveFile, int opt_iter, double obj) {
    objectiveFile << opt_iter << "," << obj << "\n";
}

void WriteConstraints(std::ofstream &constraintFile, int conLen, int opt_iter, double *result) {
    for(int j = 0; j < conLen+1; ++j)
    {
        if (j < 1) {
            constraintFile << opt_iter << ",";
        }
        else {
            constraintFile << result[j-1];
            if(j != conLen) {
                constraintFile << ","; // No comma at end of line
            }
        }
        
    }
    constraintFile << "\n";
}