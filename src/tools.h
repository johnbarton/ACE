#ifndef TOOLS_H
#define TOOLS_H


#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <math.h>
#include <stdio.h>


// Typedefs
typedef std::vector<unsigned long> Key;
typedef std::vector<std::vector<double> > Vector;
typedef std::vector<std::vector<int> > IntVector;


double computeGamma_L0(const Vector &, double);
double computeGamma_L2(const Vector &, double);

double innerProduct(double[], double[], size_t);
double innerProduct(const std::vector<double> &, const std::vector<double> &);
double innerProduct(const Vector &, const Vector &);
double innerProduct(const Vector &, const Vector &, const IntVector &);
double L1(double[], size_t);
double L1(const std::vector<double> &);
double L1(const Vector &);
double L2(double[], size_t);
double L2(const std::vector<double> &);
double L2(const Vector &);
double L2(const Vector &, const IntVector &);
double LInfinity(double[], size_t);
double LInfinity(const std::vector<double> &);
double LInfinity(const Vector &);
double LInfinity(const Vector &, const IntVector &);

int hamming(const std::vector<char> &, const std::vector<char> &);
int hamming(const std::vector<int> &, const std::vector<int> &);
bool hammingCut(const std::vector<char> &, const std::vector<char> &, int);
bool hammingCut(const std::vector<int> &, const std::vector<int> &, int);

void insertInPlace(std::vector<int> &, int);
void insertInPlace(int[], int, int);
void eraseInPlace(std::vector<int> &, int);
void eraseInPlace(int[], int, int);
void quickSort(double[][2], int, int);
void quickSort(double[][3], int, int);
void quickSort(std::vector<int> &, int, int);
bool checkSize(const Key &, const Key &);
void getSuperset(const Key &, const Key &, Key &);

void expandCorrelations(const Vector &, Vector &);
void expandCorrelations3(Vector &, Vector &);

void unpack(Vector &, std::vector<double> &, const Vector &);
void pack(Vector &, std::vector<double> &, const Vector &);
void unpack_sparse(Vector &, std::vector<double> &, const Vector &, const IntVector &);
void pack_sparse(Vector &, std::vector<double> &, const Vector &, const IntVector &);

void modifiedCholeskyInverse(double[], size_t);
void modifiedCholeskyInverse(Vector &, std::vector<double> &, int, const Vector &);
void modifiedCholeskyInverse(std::vector<double> &, std::vector<double> &, std::vector<double> &, std::vector<double> &, size_t);
void modifiedCholeskyInverse_sparse(Vector &, std::vector<double> &, int, const Vector &, const IntVector &);
void modifiedCholeskyInverse_sparse(Vector &, std::vector<double> &, std::vector<double> &, std::vector<double> &, std::vector<double> &, const Vector &, const IntVector &);
void modifiedCholeskyInverse(std::vector<double> &, int);
double modifiedCholeskyInverseDet(double[], size_t);
double modifiedCholeskyInverseDet(Vector &, int, const Vector &);
double symmetricDet(double[], int);

void givens(double, double, double &, double &);

void householder(double[], int, double[], double &);
void tridiagonalize(double[], int, std::vector<std::vector<double> > &);
void tridiagonalize(std::vector<double> &, int, std::vector<std::vector<double> > &);
void symmetricQR(double[], int, double[], std::vector<std::vector<double> > &);
void symmetricQR(std::vector<double> &, int, double[], std::vector<std::vector<double> > &);





// STRING MANIPULATION

// Converts generic to string

template <class T>

inline std::string tostring (const T& t) {
    
    std::stringstream ss;
    ss << t;
    return ss.str();
    
}


// Converts a string to an integer

inline int strtoint(const std::string &s) {
    
    std::istringstream i(s);
    int x;
    
    if (!(i >> x)) printf("String to integer conversion failed!");
    return x;
    
}


// Converts a string to a double

inline double strtodouble(const std::string &s) {
    
    std::istringstream i(s);
    double x;
    
    if (!(i >> x)) printf("String to double conversion failed!");
    return x;

}


// MISC FUNCTIONS

// Given the size of a coupling or correlation vector, return the corresponding number of spins

inline int sizetolength(size_t size) {
    
    return (int) ((sqrt(1 + 8 * size) - 1) / 2);
    
}


// Checks to see whether the ith bit of an integer in a key is occupied

inline int checkBit(unsigned long &key, int i) {
    
    return (int) (key & 1<<i);
    
}


// Flip the ith bit of an integer

inline void flipBit(unsigned long &key, int i) {
    
    key ^= 1<<i;
    
}


// For j>i, the pair {i,j} in the list {{0,1}, {0,2}, ... } is located at offset(i,length)+j

inline int offset(int i, size_t length) {
    
    return (int) (length + i * (length - 2) - (i * (i - 1)) / 2 - 1);
    
}


// For j>i, the pair {i,j} in the list {{0,0}, {0,1}, ... } is located at hindex(i,j,length)
// This routine is necessary when the diagonal term is also included

inline int hindex(int i, int j, size_t length) {
    
    return (int) (i * length - (i * (i + 1)) / 2 + j);
    
}


// Returns as above, but with a check to make sure that j>i

inline int safeHindex(int i, int j, size_t length) {
    
    if (j>i) return (int) (i * length - (i * (i + 1)) / 2 + j);
    else     return (int) (j * length - (j * (j + 1)) / 2 + i);
    
}


// Return the location of a pair of states in canonical order, {{0,0}, {0,1}, ..., {0,qj}, {1,0}, ...}

inline int sindex(int i, int j, size_t qi, size_t qj) {
    
    return (int) (i * qj + j);
    
}


// Return the location of a triple of states in canonical order

inline int sindex3(int i, int j, int k, size_t qi, size_t qj, size_t qk) {
    
    return (int) (i * qj * qk + j * qk + k);
    
}


// For j>i, the pair {i,j} in the list {{0},{1},...,{length-1},{0,1}, {0,2}, ... } is located at index(i,j,length)

inline int index(int i, int j, size_t length) {
    
    return (int)(length + i * (length - 2) - (i * (i - 1)) / 2 - 1 + j);
    
}


// For k>j>i, return the index of {i,j,k} in the list {{0,1,2},{0,1,3},...,{0,1,length-1},...}

inline int index(int i, int j, int k, int length) {
    
    return ( ( (i * (i + 1) * ( i + 2 )) / 3 
                  - (length-1) * (i * (i + 2) - 2 * j + 1) 
                  + ((int) pow(length-1, 2)) * (i - 1)
                  - (j * (j + 1) - 2 * k)                ) / 2 + (length * (length-3) / 2) );
    
}


// For l>k>j>i, return the index of {i,j,k,l} in the full list

inline int index(int i, int j, int k, int l, int length) {
    
    return ( (-i * (i + 1) * (i + 2) * (i + 3) 
              + 4 * (j * (1 + j) * (j + 2) - 3 * k * (k + 1) + 6 * l - 5 * (length-1)) 
              + 2 * (i * (11 + i * (9 + 2 * i)) - 6 * j * (j + 2) + 12 * k) * (length-1) 
              - 6 * (i * (i + 3) - 2 * j) * ((int) pow(length-1, 2))
              + 4 * (i - 1) * ((int) pow(length-1,3))                                   ) / 24 );
    
}


// Returns the sign of an integer x (and 0 if x==0)

inline int sign(int x) {
                
    static const int lookup_table[] = { -1, 1, 0, 0 };
    return lookup_table[(x==0) * 2 + (x>0)];

}


// Returns the sign of a double x (and 0 if x==0)

inline double sign(double x) {
    
    //double eps=pow(10.0,-20.0);
    
    static const double lookup_table[] = { -1, 1, 0, 0 };
    return lookup_table[(x==0) * 2 + (x>0)];
    
}


// Given components (a, b), this function computes parameters for a Givens rotation matrix G = ((c, s), (-s, c)),
// such that G^T (a, b) = (r, 0). See Matrix Computations (Golub) for details.

inline void givens(double a, double b, double &c, double &s) {
 
    if (b==0) { c=1; s=0; }
    
    else {
        
        if (fabs(b) > fabs(a)) { s = 1 / sqrt(1 + pow(a,2) / pow(b,2)); c = -a * s / b; }
        else                   { c = 1 / sqrt(1 + pow(b,2) / pow(a,2)); s = -b * c / a; }
        
    }
    
}


#endif
