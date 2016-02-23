#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "algorithm.h"  // Algorithm declarations


/*********************************************************************
 
                    COMMAND LINE INPUT FORMAT
 
 Command line instructions tell the program where to look for input
 files and where to send output, as well as the setting of various
 parameters (gamma, theta, etc) and flags (useSparse, etc).
 Note that numerical parameters may be entered in either scientific
 (recommended) or standard decimal notation. True/false switches are
 off by default - if entered in the command line, the corresponding
 option is set to true. Defaults can be found in dataStructures.h.
 The conventions are given below:
 
 -(flag name): (type of input expected)
 
 -d: string
    Default: "." (current directory)
    Path to the directory where the data file is located, and where
    output will be written.
 
 -i: string
    Default: "input"
    The location of the file containing a set of correlations
    from which to infer Ising model parameters.
 
 -o: string
    Default: "output"
    The location of the file where output is to be sent. Each
    different type of output file will have a different file type, 
    e.g. .j for couplings.
    
 -cmap: string
    Default: none 
    When the network of interactions (e.g. contact map) is known, a list
    of perselected 2-site clusters can be given. "string" represent the name 
    and the location of the file from which clusters are read. The extension 
    of the file has to be .cl
    
 -inputcl: string
    Default: none 
    When a list of interesting clusters (e.g. from previously runnings) is known,
    this list of perselected n-site clusters can be given and used for inference.
    "string" represent the name and the location of the file from which clusters are 
    read. The extension of the file has to be .cl
 
 -cl: none
    Print the list of selected clusters in a file .cl in the output folder

 -b: real number
    Default: 1.0e+4
    Number of samples used to compute the correlations. Used to 
    determine the inference error.
 
 -kmin: integer
    Default: 1
    Minimum cluster size, useful for avoiding the inference of models
    that are too sparse. The algorithm will continue to lower the
    threshold until clusters of at least this size are selected.
 
 -kmax: integer
    Default: N (system size)
    Maximum cluster size. The algorithm will halt when clusters of this
    size are selected.
 
 -t: real number
    Default: none
    Run the algorithm at the input value of theta, in scientific or
    standard decimal notation. This line is intended to be used when
    inferrence is to be done only for a single value of theta, and
    will be overridden if thetaMax and thetaMin are set different from t.
 
 -tmin: real number
    Default: 1.0e-10
    The minimum value of theta. See description of -ts below for more 
    information.
 
 -tmax: real number
    Default: 1.0e+0
    The maximum value of theta. See description of -ts below for more
    information.
 
 -ts: real number
    Default: 1.05
    The logarithmic step size to use for successive updates of theta. 
    When the program loops over different values of theta, it begins
    by running the algorithm at the largest value of the cutoff and stores 
    the cluster information. The algorithm is then re-run for successively
    smaller values of the cutoff, theta_(i+1) = theta_i / thetastep,
    until theta < thetaMin. These re-runs use the previously stored
    cluster information, so they take considerably less time to run.

 -trec: real number
    Default: 0 (any intermediate recording)
    The logarithmic step size for theta to record the inferred parameters. 
    Given this interval the chosen value corresponds the theta producing the 
    minimum error on correlations.
    
 -mcb: integer
    Default: 4.0e+4
    Number of Monte Carlo samples to take to check inference error.
    
 -mcr: integer
    Default: 1
    Number of independent Monte Carlo runs to perform.
 
 -g0: real number
    Default: 1.0e-4
    The L0 regularization strength. Using this flag also turns on L0 regularization.
 
 -g2: real number
    Default: 0.0
    The L2 regularization strength. L2 regularization is enabled by setting the 
    regularization strength to a nonzero value using this flag, or by using the
    -ag flag below.
    
 -gi: none
    Use gauge invariant L2 regularization for couplings.
 
 -ag: none
    Attempt to set the L0 and L2 regularization strengths to their optimal values,
    based on the number of samples (input) in the data. The integer input here
    *can be different* from the value used for determining the correlations.
 
 -l0: none
    If selected, L0-norm (sparse) regularization is used.
    
 -lax: none
    If selected, use a laxer cluster construction rule.

 -v: none
    Enable verbose output.
 
 *********************************************************************/


// MAIN PROGRAM

int main(int argc, char *argv[]) {
    
    RunParameters r;
    
    // Process command line input
    
    for (int i=1;i<argc;i++) {
        
        // Location of input/output files

        if      (strcmp(argv[i],"-d")==0)       { if (++i==argc) break; else r.directory=argv[i];               }
        else if (strcmp(argv[i],"-i")==0)       { if (++i==argc) break; else r.infile=argv[i];                  }
        else if (strcmp(argv[i],"-o")==0)       { if (++i==argc) break; else r.outfile=argv[i];                 }
        else if (strcmp(argv[i],"-ss")==0)      { if (++i==argc) break; else { r.useCmap=true;
                                                                               r.ssinfile=argv[i];            } }
        else if (strcmp(argv[i],"-inputcl")==0) { if (++i==argc) break; else { r.inputClusters=true;
                                                                               r.ssinfile=argv[i];            } }
        else if (strcmp(argv[i],"-cl")==0)      { r.recClusters=true;                                           }
        
        // Cluster size and threshold cutoffs
        
        else if (strcmp(argv[i],"-kmin")==0)    { if (++i==argc) break; else r.kmin=strtoint(argv[i]);          }
        else if (strcmp(argv[i],"-kmax")==0)    { if (++i==argc) break; else r.kmax=strtoint(argv[i]);          }
        else if (strcmp(argv[i],"-t")==0)       { if (++i==argc) break; else { r.thetaMax=strtodouble(argv[i]);
                                                                             r.thetaMin=strtodouble(argv[i]); } }
        else if (strcmp(argv[i],"-tmin")==0)    { if (++i==argc) break; else r.thetaMin=strtodouble(argv[i]);   }
        else if (strcmp(argv[i],"-tmax")==0)    { if (++i==argc) break; else r.thetaMax=strtodouble(argv[i]);   }
        else if (strcmp(argv[i],"-ts")==0)      { if (++i==argc) break; else r.thetaStep=strtodouble(argv[i]);  }
        else if (strcmp(argv[i],"-trec")==0)    { if (++i==argc) break; else r.recordStep=strtodouble(argv[i]); }
        
        else if (strcmp(argv[i],"-lax")==0)     { r.useLax=true;                                                }
    
        // Regularization strengths and settings
        
        else if (strcmp(argv[i],"-r")==0)       { r.useRef=true;                                                }
        else if (strcmp(argv[i],"-ag")==0)      { r.computeGamma=true;                                          }
        else if (strcmp(argv[i],"-l0")==0)      { r.useSparse=true;                                             }
        else if (strcmp(argv[i],"-g0")==0)      { if (++i==argc) break; else { r.gamma0=strtodouble(argv[i]);
                                                                               r.useSparse=true;              } }
        else if (strcmp(argv[i],"-g2")==0)      { if (++i==argc) break; else r.gamma2=strtodouble(argv[i]);     }
        else if (strcmp(argv[i],"-gi")==0)      { r.useGI=true;                                                 }
        
        // Monte Carlo settings
        
        else if (strcmp(argv[i],"-b")==0)       { if (++i==argc) break; else r.sampleB=strtodouble(argv[i]);    }
        else if (strcmp(argv[i],"-mcb")==0)     { if (++i==argc) break; else r.b=strtoint(argv[i]);             }
        else if (strcmp(argv[i],"-mcr")==0)     { if (++i==argc) break; else r.runs=strtoint(argv[i]);          } //#FLAG NOT YET IMPLEMENTED
        
        // Optional output
        
        else if (strcmp(argv[i],"-v")==0)       { r.useVerbose=true;                                            }
        
        else printf("Unrecognized command! '%s'\n",argv[i]);
                
    }
    
    run(r);
    
    return 0;
    
}