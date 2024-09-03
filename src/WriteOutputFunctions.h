#pragma once
#ifndef WRITEOUTPUTFUNCTION_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define WRITEOUTPUTFUNCTION_H

//#include "FeaHelperFunctions.h"
//#include "ProjectionHelperFunctions.h"

void WriteObjective(std::ofstream &objectiveFile, int opt_iter, double obj);
void WriteConstraints(std::ofstream &constraintFile, int conLen, int opt_iter, double *result);

#endif