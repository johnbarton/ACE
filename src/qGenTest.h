#ifndef QGENTEST_H
#define QGENTEST_H


#include <iostream>
#include <string>

#include "tools.h"  // Numerical tools


// This class holds the parameters needed for running the model validation routine

class RunParameters {
    
public:
    
    std::string directory;  // Path to the directory where the inut file is located
                            // Output will also be sent to this directory
    std::string infile;     // Input file from which couplings are to be read
    std::string consfile;   // Input file from which reference sequence is to be read
    std::string weightfile; // Input file from which weight vector is to be read
    std::string msafile;    // Input MSA file 
    std::string outfile;    // Output file (prefix) where data is to be written
    
    double b;               // Number of data points to take in a single run
    double runs;            // Number of times to run dynamics
    
    double sampleB;         // Number of data points in the sample
    bool p3;                // compute 3-point correlations
    bool p3red;             // compute 3-point correlations and only some small 3-point correlations will be printed
    double gamma;           // Regularization strength
    bool useGamma;          // If true, use L2 regularization (with strength computed using sampleB and correlations)
    bool useGI;             // If true, use gauge invariant regularization for couplings
    
    bool useStart;          // Use different starting sequence than all wild-type
    bool refMC;             // Use reference Monte Carlo
    std::string startfile;  // File containing the starting sequence
    bool ThreePoints;       // compute 3 point correlations
    bool MSAout;            // printf MonteCarlo alignment and compute energies
    
    bool useVerbose;        // If true, print extra information while program is running
    
    RunParameters() {
        
        directory=".";
        infile="input";
        consfile="input";
        weightfile="input";
        msafile="input";
        outfile="output";
        b=800000;
        runs=1;
        
        sampleB=1000;
        p3=false;
        p3red=false;
        gamma=0;
        useGamma=false;
        useGI=false;
        
        useStart=false;
        refMC=false;
        useVerbose=false;
    
    }
    
    std::string getInfile()          { return (directory+"/"+infile+".j");       }
    std::string getStartInfile()     { return (directory+"/"+startfile);         }
    std::string getInfileAl()        { return (directory+"/"+msafile+".cmsa");   }
    std::string getConsensusInfile() { return (directory+"/"+consfile+".cons");  }
    std::string getWeights()         { return (directory+"/"+weightfile+".wgt"); }
    std::string getMOutfile() 	     { return (directory+"/"+outfile+".m");      }
    std::string getP2Outfile()       { return (directory+"/"+outfile+".p2");     }
    std::string getCCOutfile()       { return (directory+"/"+outfile+".c2");     }
    std::string getPKOutfile()       { return (directory+"/"+outfile+".pk");     }
    std::string getP3Outfile()       { return (directory+"/"+outfile+".p3");     }
    std::string getC3Outfile()       { return (directory+"/"+outfile+".c3");     }
    
    ~RunParameters() {}
    
};

void runGenTest(RunParameters &);


#endif
