#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H


#include <vector>
#include <math.h>
#include <string>
#include <sstream>
#include <cstring>
#include <iostream>
#include <bitset>
#include <stdio.h>

#include "tools.h"  // Numerical tools


// CLUSTER CLASSES

// Cluster class, which holds the information for a given cluster

class Cluster {
    
public:
	
	double dS;
	Vector dJ;
    bool selected;
    bool superSelected;
	
	Cluster() {}
    Cluster(int size) {
        
		dS = 0;
        dJ.resize((size * (size + 1)) / 2);
        selected = false;
        superSelected = false;
        
    }
    ~Cluster() {}
	
};


// Functor which is used to determine whether or not a cluster should
// be counted among the significant clusters

class Significant {
    
public:
    
    double theta;
    
    Significant(double t) : theta(t) {}
    
    bool operator() (const Cluster &c) {
        
        return ((fabs(c.dS) > theta) || (c.dJ.size()==1) || c.selected);
        
    }
    
};


// PROGRAM SETTINGS

// This class holds the parameters needed for running the adaptive cluster algorithm

class RunParameters {
    
public:
    
    std::string directory;  // Path to the directory where the inut file is located
                            // Output will also be sent to this directory
    std::string infile;     // Input file (prefix)
    std::string ssinfile;   // Input file for secondary structure (default selected clusters)
    std::string outfile;    // Output file (prefix)
    
    int kmin;               // Minimum cluster size before computation is truncated
    int kmax;               // Maximum cluster size before computation is truncated
    int knmaxk;             // k parameter for cluster truncation rule (stop when >= n clusters of size k found)
    int knmaxn;             // n parameter for cluster truncation rule (stop when >= n clusters of size k found)
    double sampleB;         // Number of samples used to compute correlations from data
    double gamma0;          // Reference entropy regularization strength
    double gamma2;          // Cluster entropy regularization strength
    double gammah;          // L2 regularization strength multiplier for fields
    double thetaMin;        // The minimum cutoff value (loop starts here)
    double thetaMax;        // The maximum cutoff value (loop ends here)
    double recordStep;      // Step size when recording final data
    double thetaStep;       // Size of logarithmic steps when looping over theta
    int b;                  // Number of data points to take in a single Monte Carlo run
    int runs;               // Number of times to run Monte Carlo dynamics
    bool useGI;             // If true, use gauge invariant regularization for couplings
    bool computeGamma;      // If true, compute regularization strength based on number of samples
    bool useLax;            // If true, use a laxer cluster construction rule
    bool useRef;            // If false, do not use the reference entropy S0 (i.e. set S0, J0 to zero)
    bool useSparse;         // If true, use sparse L0-norm regularization, in addition to the default L2
    bool useCmap;           // If true, the inference is performed around the give contact map
    bool inputClusters;     // If true, select clusters from a given list and perform inference
    bool recClusters;       // If true, save the list of cluster selected at convergence in a .wuss file
    bool recClusterCover;   // If true, record all selected clusters that are not subsets of
                            // other selected clusters in a file ending in -co.dat
    bool useVerbose;        // If true, print extra information while program is running
    bool useVeryVerbose;    // If true, use very verbose output

    
    RunParameters() {
        
        directory = ".";
        infile    = "input";
        ssinfile  = "input";
        outfile   = "output";
        
        kmin       = 0;
        kmax       = 0;
        knmaxk     = 0;
        knmaxn     = 1;
        sampleB    = 1000;
        gamma0     = 1.0e-4;
        gamma2     = 0.0;
        gammah     = 0.01;
        thetaMin   = 1.0e-10;
        thetaMax   = 1.0e+0;
        thetaStep  = 1.05;
        recordStep = 0;
        b          = 40000;
        runs       = 1;
        
        useGI           = false;
        computeGamma    = false;
        useLax          = false;
        useRef          = false;
        useSparse       = false;
        useCmap         = false;
        inputClusters   = false;
        recClusters     = false;
        recClusterCover = false;
        useVerbose      = false;
        useVeryVerbose  = false;
        
    }
    std::string getCorrelationsInfile()   { return (directory+"/"+infile+".p");       }
    std::string getSecStructInfile()      { return (directory+"/"+ssinfile+".cl");    }
    std::string getClusterOutfile()       { return (directory+"/"+outfile+".cl");     }
    std::string getClusterCoverOutfile()  { return (directory+"/"+outfile+"-cc.dat"); }
    std::string getCouplingsOutfile()     { return (directory+"/"+outfile+".j");      }
    std::string getBestCouplingsOutfile() { return (directory+"/"+outfile+"-best.j"); }
    std::string getOutfile_TH()           { return (directory+"/"+outfile);           }  // Theta recording files
    std::string getCorrelationsOutfile()  { return (directory+"/"+outfile+".p");      }  //#FLAG
    std::string getSupplementaryOutfile() { return (directory+"/"+outfile+".sce");    }
    ~RunParameters() {}
    
};


#endif
