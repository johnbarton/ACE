#ifndef SPARSEINVERSE_H
#define SPARSEINVERSE_H


#include "inverse.h"    // Regularization functions


// Partition function
double getZ(int, const Vector &, int[], int[], int[], Vector &, Vector &, int, double, double[], const IntVector &, const IntVector &);
void optimizeS(const Vector &, const Vector &, double &, Vector &, Vector &, const Vector &, const IntVector &, const IntVector &);
double getZ_gradOnly(int, const Vector &, int[], int[], int[], Vector &, int, double, double[], const IntVector &, const IntVector &);
void optimizeS_gradOnly(const Vector &, const Vector &, double &, Vector &, const Vector &, const IntVector &, const IntVector &);

// Auxiliary Functions
void computeS0andJ0_Empty(const Vector &, int, double, Vector &, double &);
void initialize(const Vector &, int, Vector &, Vector &, int &, const IntVector &, IntVector &, IntVector &);
bool useDescent(const Vector &, const IntVector &);
void computeDescentStep(const Vector &, int, Vector &, const IntVector &);
void computeNewtonStep(const Vector &, Vector &, std::vector<double> &, std::vector<double> &, std::vector<double> &, std::vector<double> &, Vector &, const IntVector &);
double lineSearch_simple(const Vector &, const Vector &, Vector &, double, double, const Vector &, Vector &, Vector &, Vector &, double, bool, double, const IntVector &, const IntVector &, const IntVector &);
double lineSearch_simple(const Vector &, const Vector &, Vector &, double, double, const Vector &, Vector &, Vector &, Vector &, double, const IntVector &, const IntVector &, const IntVector &);
double lineSearch_interp(const Vector &, const Vector &, Vector &, double, double, const Vector &, Vector &, Vector &, Vector &, double, const IntVector &, const IntVector &, const IntVector &);
void makeStep(Vector &, const Vector &, const Vector &, double, double, double, Vector &, Vector &, double &, const IntVector &, const IntVector &, const IntVector &);

//L2 regularization
void computeSandJ_L2_sparse(const Vector &, int, double, double, Vector &, double &, const IntVector &, Vector &, Vector &, Vector &, Vector &, Vector &, Vector &, std::vector<double> &, std::vector<double> &, std::vector<double> &, std::vector<double> &);

// L0 regularization
void computeS0andJ0_L0(const Vector &, int, double, Vector &, double &);
void computeSandJ_L0(const Vector &, int, double, double, double, Vector &, double &);
void regularizeS_L2(const Vector &, double &, Vector &, Vector &, const Vector &, double, double, const IntVector &);
void regularizeS_L2_gradOnly(const Vector &, double &, Vector &, const Vector &, double, double, const IntVector &);


#endif
