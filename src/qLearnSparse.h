#ifndef QLEARNSPARSE_H
#define QLEARNSPARSE_H


#include <iostream>
#include <string>

#include "tools.h"  // Numerical tools


// This class holds the parameters needed for running the ising program

class RunParameters {
    
public:
    
    std::string directory;  // Path to the directory where the inut file is located
                            // Output will also be sent to this directory
    std::string infile;     // Input file from which couplings are to be read
    std::string compfile;   // Input file from which correlations are to be read for comparisons
    std::string outfile;    // Output file (prefix) where data is to be written
    
    double b;               // Number of data points to take in a single run
    double runs;            // Number of times to run dynamics
    double epsilon;         // Acceptable error between target and current correlations
    
    double sampleB;         // Number of data points in the sample
    double gamma;           // Regularization strength
    bool useGamma;          // If true, use L2 regularization (with strength computed using sampleB and correlations)
    bool useGI;             // If true, use gauge invariant regularization for couplings

    double a;               // Line search acceleration multiplier
    double d;               // Line search deceleration multiplier
    double s;               // Line search default step size
    bool useSparse;         // If true, maintain sparse couplings (do not add new couplings)
    
    bool useStart;          // Use different starting sequence than all wild-type
    std::string startfile;  // File containing the starting sequence
    
    bool useVerbose;        // If true, print extra information while program is running
    
    RunParameters() {
        
        directory=".";
        infile="input";
        outfile="output";
        compfile="input";
        b=800000;
        runs=1;
        epsilon=1.0;
        
        sampleB=1000;
        gamma=0;
        useGamma=false;
        useGI=false;
        
        a=2.0;
        d=0.5;
        s=0.001;
        useSparse=true;
        
        useStart=false;
        useVerbose=false;
        
    }
    std::string getInfile()              { return (directory+"/"+infile+".j");      }
    std::string getCompareInfile()       { return (directory+"/"+compfile+".p");    }
    std::string getCompareOutfile()      { return (directory+"/"+outfile+".fit");   }
    std::string getCorrelationsOutfile() { return (directory+"/"+outfile+".p");     }
    std::string getCouplingsOutfile()    { return (directory+"/"+outfile+".j");     }
    std::string getStartInfile()         { return (directory+"/"+startfile+".dat"); }
    ~RunParameters() {}
    
};


// Monte Carlo
void updateStep(   const Vector &, const Vector &, const IntVector &, double, double, Vector &, Vector &, Vector &, Vector &);
void updateStep_GI(const Vector &, const Vector &, const IntVector &, double, double, Vector &, Vector &, Vector &, Vector &);
bool isClose(double, double, double);
bool isCloseC(double, double, double, double, double);
void chop(Vector &);
void runLearn(RunParameters &);


#endif
