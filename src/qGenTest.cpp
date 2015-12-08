#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <cstring>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>

#include "qGenTest.h"       // Generative Test declarations
#include "monteCarlo.h"     // Monte Carlo
#include "tools.h"          // Numerical tools
#include "io.h"             // Data I/O
#include "mtrnd.h"          // Random number generator


// Runs the main program

void runGenTest(RunParameters &r) {

    // Define variables
    
    Vector J, expJ;
    std::vector<int> cons;
    std::vector<double> weight;
    
    if (r.useGI) {
    
        epsilonP2_ptr=&epsilonP2_GI;
        epsilonC_ptr=&epsilonC_GI;
        getMaxError_ptr=&getMaxError_GI;
        
    }
    else {
    
        epsilonP2_ptr=&epsilonP2;
        epsilonC_ptr=&epsilonC;
        getMaxError_ptr=&getMaxError;
        
    }
    
    // Get reference sequence from file
    
    FILE *consIn = fopen(r.getConsensusInfile().c_str(),"r");
         
    if (consIn!=NULL) getConsensus(consIn,cons);
    else { printf("Error reading input from file %s\n\n",r.getConsensusInfile().c_str()); exit(1); }
    fclose(consIn);
    
    if (r.useVerbose) {
      
        printf("Reference sequence: ");
        for (int i=0;i<cons.size();i++) printf(" %d",cons[i]);
        printf("\n\n");
	
    }

    // Retrieve couplings from file
    
    FILE *dataIn=fopen(r.getInfile().c_str(),"r");
        
    if (dataIn!=NULL) getCouplings(dataIn,J);
    else { printf("Error reading input from file %s",r.getInfile().c_str()); exit(1); }
    fclose(dataIn);

    // Resize expJ
    
    for (int i=0;i<J.size();i++) expJ.push_back(std::vector<double>(J[i].size(),0));
    for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) expJ[i][j] = exp(J[i][j]); }
    
    // Declare 2-point correlations, 3-point correlations, P(k) and magnetisations
    
    bool ThreePoints = (r.p3red || r.p3);
    int N = sizetolength(J.size()); // System size
    double alpha = 0.01;            // Field regularization multiplier
    double gamma = 0;               // Regularization strength (L2, set below)
    
    if (r.useGamma) {
    
        if (r.gamma==0) gamma=1/(r.sampleB);
        else            gamma=r.gamma;
        
    }
    
    Vector p(J.size(),std::vector<double>());           // MC magnetisations and 2-point correlations
    Vector cc(J.size(),std::vector<double>());          // MC connected 2-point correlations
    Vector q(J.size(),std::vector<double>());           // MSA magnetisations and 2-point correlations
    Vector qcc(J.size(),std::vector<double>());         // MSA connected 2-point correlations
    std::vector<std::vector<std::vector<std::vector<double> > > > p3(N);  // MC 3-point correlations
    std::vector<std::vector<std::vector<std::vector<double> > > > c3(N);  // MC connected 3-point correlations
    std::vector<std::vector<std::vector<std::vector<double> > > > q3(N);  // MSA 3-point correlations
    std::vector<std::vector<std::vector<std::vector<double> > > > qc3(N); // MSA connected 3-point correlations
    std::vector<double> pk(N+1,0);                      // MC mutation probability
    std::vector<double> qk(N+1,0);                      // MSA mutation probability    
    std::vector<double> absErr(2,0);                    // Absolute errors on magnetisation and 2-point correlations
    
    for (int i=0;i<J.size();i++) {
    
       cc[i].resize(J[i].size(),0); 
       p[i].resize(J[i].size(),0);
       qcc[i].resize(J[i].size(),0); 
       q[i].resize(J[i].size(),0);
       
    }
	if (ThreePoints) { for (int i=0;i<N;i++) {
    
	    p3[i].resize(N);
	    c3[i].resize(N);
	    q3[i].resize(N);
	    qc3[i].resize(N);
        
	    for (int j=0;j<N;j++) {
        
            p3[i][j].resize(N);
            c3[i][j].resize(N);
            q3[i][j].resize(N);
            qc3[i][j].resize(N);
            
            for (int k=0;k<N;k++) {
            
                p3[i][j][k].resize(p[i].size()*p[j].size()*p[k].size(),0);
                c3[i][j][k].resize(p[i].size()*p[j].size()*p[k].size(),0);
                q3[i][j][k].resize(p[i].size()*p[j].size()*p[k].size(),0);
                qc3[i][j][k].resize(p[i].size()*p[j].size()*p[k].size(),0);
                
            }
            
	    }
        
	} }

    // Get sequences from MSA file and compute correlations
    
    FILE *alIn=fopen(r.getInfileAl().c_str(),"r");
    FILE *weightIn=fopen(r.getWeights().c_str(),"r");

    if (alIn!=NULL){
        if (ThreePoints) getAlignment(alIn,weightIn,J,q,q3,qk,cons);
        else getAlignment(alIn,weightIn,J,q,qk,cons);
    }
    else { printf("Error reading input from file %s\n\n",r.getInfileAl().c_str()); exit(1); }
    fclose(alIn);    
    if (weightIn!=NULL) fclose(weightIn);
        
    if (r.useVerbose) printf("Got N=%d, len(h[0])=%d\n",N,(int)J[0].size());
        
    // Get default starting configuration, if nontrivial
    
    std::vector<int> lattice(N);
    
    if (r.useStart) {
    
        FILE *startIn=fopen(r.getStartInfile().c_str(),"r");
        for (int i=0;i<N;i++) fscanf(startIn,"%d",&lattice[i]);
    
    }
    
    else { for (int i=0;i<N;i++) lattice[i]=(int) p[i].size(); } 
    
    // Prepare to simulate
    
    srand((unsigned)time(0));

    // Run MC and get correlations
   
    if (ThreePoints) getErrorGenTest(J, expJ, r.sampleB, r.b, r.runs, p, lattice, pk, p3, cons); // compute errors on P P2 and MAX
    else             getErrorGenTest(J, expJ, r.sampleB, r.b, r.runs, p, lattice, pk, cons);     // compute errors on P P2 and MAX

    //Compute connected correlations
    
    double Neff  = 0;
    double NJeff = 0;
    // estimate the threshold for correlations to print out
    double meanq = 0;
    
    for (int i=0;i<lattice.size();i++) { for (int a=0;a<p[i].size();a++) {

        Neff++;
        meanq+=q[i][a];

        absErr[0] += (p[i][a] - q[i][a]) * (p[i][a] - q[i][a]);
    
        for (int j=i+1;j<lattice.size();j++) { for (int b=0;b<p[j].size();b++) {

            NJeff++;
	    
            int idx = index(i,j,lattice.size());
            int sab = sindex(a,b,J[i].size(),J[j].size());
            
            absErr[1] += (p[idx][sab] - q[idx][sab]) * (p[idx][sab] - q[idx][sab]);

            cc[idx][sab]  = p[idx][sab] - (p[i][a] * p[j][b]);
            qcc[idx][sab] = q[idx][sab] - (q[i][a] * q[j][b]);
		    
            if (ThreePoints) { for (int k=j+1;k<lattice.size();k++) { for (int c=0;c<p[k].size();c++) {
		
                int ijx  = idx;
                int ikx  = index(i,k,lattice.size());
                int jkx  = index(j,k,lattice.size());
                int sac  = sindex(a,c,J[i].size(),J[k].size());
                int sbc  = sindex(b,c,J[j].size(),J[k].size());
                int sabc = sindex3(a,b,c,J[i].size(),J[j].size(),J[k].size());
                
                c3[i][j][k][sabc]  = p3[i][j][k][sabc] - (p[i][a]*p[jkx][sbc]) - (p[j][b]*p[ikx][sac]) - (p[k][c]*p[ijx][sab]) + (2*(p[i][a]*p[j][b]*p[k][c]));
                qc3[i][j][k][sabc] = q3[i][j][k][sabc] - (q[i][a]*q[jkx][sbc]) - (q[j][b]*q[ikx][sac]) - (q[k][c]*q[ijx][sab]) + (2*(q[i][a]*q[j][b]*q[k][c]));
		    
            } } }
            
        } }
        
    } }
	
    absErr[0] = sqrt(absErr[0]/Neff);
    absErr[1] = sqrt(absErr[1]/NJeff);
    meanq=meanq/Neff;
    
    // Print out errors

    double maxPrecision=1/(r.sampleB);

    double ep1 = epsilonP(q, p, N, maxPrecision, J, gamma, alpha);
    double ep2 = (*epsilonP2_ptr)(q, p, N, maxPrecision, J, gamma);
    double em  = (*getMaxError_ptr)( q, p, maxPrecision, J, gamma, alpha);
    
    printf("\nRelative errors: P %f, P2 %f MAX %f gamma %f\n",ep1,ep2,em,gamma);
    printf("Absolute errors: P %f, P2 %f \n\n",absErr[0],absErr[1]);

    //Print results for comparison
    
    FILE *mOut  = fopen(r.getMOutfile().c_str(),"w");
    FILE *pOut  = fopen(r.getP2Outfile().c_str(),"w");
    FILE *ccOut = fopen(r.getCCOutfile().c_str(),"w");
    FILE *pkOut = fopen(r.getPKOutfile().c_str(),"w");
    
    printMagnetisations(mOut, q, p);
    double num=0;
    printCorrelations(ccOut, qcc, cc,pOut, q, p);
    
    if (ThreePoints){
      
        FILE *p3Out = fopen(r.getP3Outfile().c_str(),"w");
        FILE *c3Out = fopen(r.getC3Outfile().c_str(),"w");
        if (r.p3red) num=meanq*meanq*meanq;
        if (r.p3) num=0;
        print3points(c3Out, qc3, c3,p3Out, q3, p3, num);
	
    }
    for (int sit=0;sit<N;sit++) fprintf(pkOut,"%d %le %le\n",sit,qk[sit],pk[sit]);
    fflush(pkOut);

}

