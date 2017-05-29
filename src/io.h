#ifndef IO_H
#define IO_H


#include <vector>
#include <math.h>
#include <stdio.h>

#include "tools.h"      // Numerical tools


void getCorrelations(FILE *, Vector &);
void getCouplings(FILE *, Vector &);
void getSS(FILE *, IntVector &p);
void getMSA(FILE *, std::vector<std::vector<char> > &, std::vector<std::string> &);
void printMSA(FILE *, const std::vector<std::vector<char> > &, const std::vector<std::string> &);
void printCouplings(FILE *, const Vector &);
void printSupplementaryOutput(FILE *, double, const std::vector<double> &, double, int, unsigned long, unsigned long);

void getConsensus(FILE *, std::vector<int> &);
void getWeights(FILE *, std::vector<double> &);
void getAlignment(FILE *, FILE *, Vector &, Vector &, Vector &, std::vector<double> &, std::vector<int> &, std::string);
void getAlignment(FILE *, FILE *, Vector &, Vector &, std::vector<double> &, std::vector<int> &, std::string);
void printMagnetisations(FILE *, const Vector &, const Vector &);
void printCorrelations(FILE *, const Vector &, const Vector &,FILE *, const Vector &, const Vector &,double);
void printCorrelationsError(FILE *, const Vector &, const Vector &, const Vector &, FILE *, const Vector &, const Vector &, const Vector &, double);
void print3points(FILE *, const std::vector<std::vector<std::vector<std::vector<double> > > > &, const std::vector<std::vector<std::vector<std::vector<double> > > > &, FILE *, const std::vector<std::vector<std::vector<std::vector<double> > > > &, const std::vector<std::vector<std::vector<std::vector<double> > > > &, double);
void print3pointserror(FILE *, const std::vector<std::vector<std::vector<std::vector<double> > > > &, const std::vector<std::vector<std::vector<std::vector<double> > > > &, const std::vector<std::vector<std::vector<std::vector<double> > > > &, FILE *, const std::vector<std::vector<std::vector<std::vector<double> > > > &, const std::vector<std::vector<std::vector<std::vector<double> > > > &, const std::vector<std::vector<std::vector<std::vector<double> > > > &, double);


#endif
