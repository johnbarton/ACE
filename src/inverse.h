#ifndef INVERSE_H
#define INVERSE_H


// GLOBAL VARIABLES

extern void (*computeS0andJ0_ptr)(const Vector &, int, double, Vector &, double &);
extern void (*computeSandJ_ptr)(const Vector &, int, double, double, double, Vector &, double &);
extern void (*regularizeS_ptr)(const Vector &, double &, Vector &, Vector &, const Vector &, double, double);
extern void (*regularizeS_gradOnly_ptr)(const Vector &, double &, Vector &, const Vector &, double, double);


// Partition function
double getZ(int, const Vector &, std::vector<int> &, std::vector<int> &, Vector &, Vector &, int, double);
double getZ_gradOnly(int, const Vector &, std::vector<int> &, std::vector<int> &, Vector &, int, double);
void optimizeS(const Vector &, const Vector &, double &, Vector &, Vector &, const Vector &);
void optimizeS_gradOnly(const Vector &, const Vector &, double &, Vector &, const Vector &);

// Auxiliary functions
void computeS0andJ0_Empty(const Vector &, int, double, Vector &, double &);
void initialize(const Vector &, int, Vector &, Vector &, int &);
bool useDescent(const Vector &);
void computeDescentStep(const Vector &, int, Vector &);
void computeNewtonStep(const Vector &, Vector &, std::vector<double> &, int, int, Vector &);
//double lineSearch_simple(const Vector &, const Vector &, Vector &, double, const Vector &, Vector &, Vector &, Vector &, double);
double lineSearch_simple(const Vector &, const Vector &, Vector &, double, double, const Vector &, Vector &, Vector &, Vector &, double, bool, double);
double lineSearch_interp(const Vector &, const Vector &, Vector &, double, double, const Vector &, Vector &, Vector &, Vector &, double);
void makeStep(Vector &, const Vector &, const Vector &, double, double, double, Vector &, Vector &, double &);

// Main functions
void computeS0andJ0_L2(const Vector &, int, double, Vector &, double &);
void computeS0andJ0_L2_binary(const Vector &, int, double, Vector &, double &);
void computeSandJ_L2(const Vector &, int, double, double, double, Vector &, double &);
void regularizeS_L2(const Vector &, double &, Vector &, Vector &, const Vector &, double, double);
void regularizeS_L2_gradOnly(const Vector &, double &, Vector &, const Vector &, double, double);
void regularizeS_L2_GI(const Vector &, double &, Vector &, Vector &, const Vector &, double, double);
void regularizeS_L2_gradOnly_GI(const Vector &, double &, Vector &, const Vector &, double, double);


#endif
