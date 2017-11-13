#ifndef ALGORITHM_H
#define ALGORITHM_H


#include <set>
#include <map>

#include "dataStructures.h" // Data structures
#include "tools.h"          // Numerical tools


// Cluster-related
void makeSingleCluster(Cluster &, const Key &, int);
void makeCluster(Cluster &, const Key &, int);
void computeDSandDJ(int[], int, Vector &, double &);
double computeGamma_L2(const Vector &, double);

// Main algorithm
void selectClusters(std::set<Key> &, int, int, int, int, double, bool, FILE *);
void getClusters(std::set<Key> &, int &, int, double, bool, bool, const std::vector<int> &, FILE *);
void getCouplings(Vector &, double &, unsigned long &, bool &, int, int, double);
int run(RunParameters &);


#endif
