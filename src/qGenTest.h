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
    double gamma;           // Regularization strength
    double nmax;            // (Approximate) maximum number of 3-point correlations to record
    double pthresh;         // Record 3-point correlations larger than this size in absolute value
    bool useGamma;          // If true, use L2 regularization (with strength computed using sampleB and correlations)
    bool useGI;             // If true, use gauge invariant regularization for couplings
    bool computeP2;         // If true, compute 2-point correlations
    bool computeP3;         // If true, compute 3-point correlations
    bool useNMax;           // If true, specify a maximum number of 3-point correlations to write to file
    bool recMSA;            // If true, print Monte Carlo alignment and compute energies
    bool useVerbose;        // If true, print extra information while program is running
    bool useStart;          // Use different starting sequence than all wild-type
    std::string startfile;  // File containing the starting sequence
    
    
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
        gamma=0;
        nmax=0;
        pthresh=1e-4;
        useGamma=false;
        useGI=false;
        computeP2=true;
        computeP3=false;
        useNMax=false;
        recMSA=false;
        useVerbose=false;
        useStart=false;
    
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
    std::string MSAOutfile()         { return (directory+"/"+outfile+".mc");     } // msa -> mc
    std::string EnergiesOutfile()    { return (directory+"/"+outfile+".mce");    } // e -> mce
    std::string MSAEnOutfile()       { return (directory+"/"+outfile+".msae");   }
    
    ~RunParameters() {}
    
};


int runGenTest(RunParameters &);


#endif
