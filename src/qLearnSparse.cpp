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

#include "qLearnSparse.h"   // MC learning declarations
#include "monteCarlo.h"     // Monte Carlo
#include "tools.h"          // Numerical tools
#include "io.h"             // Data I/O
#include "mtrnd.h"          // Random number generator


#define STEPS 1.0e-4
#define STEPA 1.9e+0
#define STEPD 5.0e-1

// GLOBAL VARIABLES

double alpha=0.0;

void (*updateStep_ptr)(const Vector &, const Vector &, const IntVector &, double, double, Vector &, Vector &, Vector &, Vector &);





// Update coupling-specific weights using a variant of the RPROP algorithm

void updateStep(const Vector &q, const Vector &p, const IntVector &nz, double gamma, double B, Vector &J, Vector &expJ, Vector &grad, Vector &weight) {

    int N            = sizetolength(J.size());
    double weightMax = 1.0 / STEPS;
    double weightMin = 0.1;
    double stepMax   = 0.1 / STEPS;
    double Jcut      = 15.0;
    double gammaexp  = gamma * Jcut * exp(-Jcut) / 2;
    
     //Update fields

    for (int i=0;i<N;i++) { for (int a=0;a<q[i].size();a++) {
        
        double K       = J[i][a];
        double newGrad = p[i][a] - q[i][a] + (2 * alpha * gamma * K);
            
        // Exponential penalty on very large fields
        if (fabs(K)>Jcut) newGrad += (alpha * gammaexp * sign(K) * exp(fabs(K))) - (2 * alpha * gamma * K);
            
        // Only perform step if correlation differs significantly from correct
        //if (isClose(q[i][a], newGrad, B)) weight[i][a] = 1;
        //
        //else {
            
            // If new gradient has same sign, accelerate, else decelerate
            if      (grad[i][a] * newGrad >= 0) weight[i][a] *= STEPA;
            else if (weight[i][a] < weightMin)  weight[i][a]  = weightMin;
            else                                weight[i][a] *= STEPD;
            
            // Sanity check
            if (weight[i][a] > weightMax)               weight[i][a] = weightMax;
            if (fabs(newGrad * weight[i][a]) > stepMax) weight[i][a] = stepMax / fabs(newGrad);
            
            J[i][a] -= newGrad * STEPS * weight[i][a];
            
        //}
        
        expJ[i][a] = exp(J[i][a]);
        grad[i][a] = newGrad;
    
    } }
    
    // Update (nonzero) couplings
    
    for (int n=0;n<nz.size();n++) {
    
        int i   = nz[n][0];
        int j   = nz[n][1];
        int idx = index(i,j,N);
        
        for (int a=0;a<J[i].size();a++) { for (int b=0;b<J[j].size();b++) {

            int    sab     = sindex(a,b,J[i].size(),J[j].size());
            double K       = J[idx][sab];
            double newGrad = p[idx][sab] - q[idx][sab] + (2 * gamma * K);
            
            // Exponential penalty on very large couplings
            if (fabs(K)>Jcut) newGrad += (gammaexp * sign(K) * exp(fabs(K))) - (2 * gamma * K);
            
            // Only perform step if correlation differs significantly from correct
            //if (isClose(q[idx][sab], newGrad, B)) weight[idx][sab] = 1;
            //
            //else {
    
                // If new gradient has same sign, accelerate, else decelerate
                if      (grad[idx][sab] * newGrad >= 0) weight[idx][sab] *= STEPA;
                else if (weight[idx][sab] < weightMin)  weight[idx][sab]  = weightMin;
                else                                    weight[idx][sab] *= STEPD;
                
                // Sanity check
                if (weight[idx][sab] > weightMax)               weight[idx][sab] = weightMax;
                if (fabs(newGrad * weight[idx][sab]) > stepMax) weight[idx][sab] = stepMax / fabs(newGrad);
            
                J[idx][sab] -= newGrad * STEPS * weight[idx][sab];
            
            //}
            
            expJ[idx][sab] = exp(J[idx][sab]);
            grad[idx][sab] = newGrad;

        } }
    
    }

}


// Update coupling-specific weights using a variant of the RPROP algorithm

void updateStep_GI(const Vector &q, const Vector &p, const IntVector &nz, double gamma, double B, Vector &J, Vector &expJ, Vector &grad, Vector &weight) {

    int N            = sizetolength(J.size());
    double weightMax = 1.0 / STEPS;
    double weightMin = 0.1;
    double stepMax   = 0.1  / STEPS;
    double Jcut      = 15.0;
    double gammaexp  = gamma * Jcut * exp(-Jcut) / 2;
    
    //Update fields

    for (int i=0;i<N;i++) { for (int a=0;a<q[i].size();a++) {
        
        double K       = J[i][a];
        double newGrad = p[i][a] - q[i][a] + (2 * alpha * gamma * K);
            
        // Exponential penalty on very large fields
        if (fabs(K)>Jcut) newGrad += (alpha * gammaexp * sign(K) * exp(fabs(K))) - (2 * alpha * gamma * K);
            
        // Only perform step if correlation differs significantly from correct
        //if (isClose(q[i][a], newGrad, B)) weight[i][a] = 1;
        //
        //else {
            
            // If new gradient has same sign, accelerate, else decelerate
            if      (grad[i][a] * newGrad >= 0) weight[i][a] *= STEPA;
            else if (weight[i][a] < weightMin)  weight[i][a]  = weightMin;
            else                                weight[i][a] *= STEPD;
            
            // Sanity check
            if (weight[i][a] > weightMax)               weight[i][a] = weightMax;
            if (fabs(newGrad * weight[i][a]) > stepMax) weight[i][a] = stepMax / fabs(newGrad);
            
            J[i][a] -= newGrad * STEPS * weight[i][a];
            
        //}
        
        expJ[i][a] = exp(J[i][a]);
        grad[i][a] = newGrad;
    
    } }
    
    // Update (nonzero) couplings
    
    for (int n=0;n<nz.size();n++) {
    
        int i   = nz[n][0];
        int j   = nz[n][1];
        int idx = index(i,j,N);
        
        // Get quantities for computing K_ij^ab

        double inv_qi = 1/((double)J[i].size()+1);
        double inv_qj = 1/((double)J[j].size()+1);
    
        double suma[J[i].size()];
        double sumb[J[j].size()];
        double sumall=0;
        
        for (int a=0;a<J[i].size();a++) suma[a]=0;
        for (int b=0;b<J[j].size();b++) sumb[b]=0;
        
        for (int a=0;a<J[i].size();a++) { for (int b=0;b<J[j].size();b++) {
            
            suma[a] += J[idx][sindex(a,b,J[i].size(),J[j].size())];
            sumb[b] += J[idx][sindex(a,b,J[i].size(),J[j].size())];
            
        } }
        
        for (int a=0;a<J[i].size();a++) sumall += suma[a];
        
        for (int a=0;a<J[i].size();a++) { for (int b=0;b<J[j].size();b++) {

            int    sab     = sindex(a,b,J[i].size(),J[j].size());
            double K       = J[idx][sab] - (suma[a] * inv_qj) - (sumb[b] * inv_qi) + (sumall * inv_qi * inv_qj);
            double newGrad = p[idx][sab] - q[idx][sab] + (2 * gamma * K);
            
            // Exponential penalty on very large couplings
            if (fabs(K)>Jcut) newGrad += (gammaexp * sign(K) * exp(fabs(K))) - (2 * gamma * K);
            
            // Only perform step if correlation differs significantly from correct
            //if (isClose(q[idx][sab], newGrad, B)) weight[idx][sab] = 1;
            //
            //else {
    
                // If new gradient has same sign, accelerate, else decelerate
                if      (grad[idx][sab] * newGrad >= 0) weight[idx][sab] *= STEPA;
                else if (weight[idx][sab] < weightMin)  weight[idx][sab]  = weightMin;
                else                                    weight[idx][sab] *= STEPD;
                
                // Sanity check
                if (weight[idx][sab] > weightMax)               weight[idx][sab] = weightMax;
                if (fabs(newGrad * weight[idx][sab]) > stepMax) weight[idx][sab] = stepMax / fabs(newGrad);
            
                J[idx][sab] -= newGrad * STEPS * weight[idx][sab];
            
            //}
            
            expJ[idx][sab] = exp(J[idx][sab]);
            grad[idx][sab] = newGrad;

        } }
    
    }

}


// Evaluate whether or not a correlation is close to target value

bool isClose(double q, double grad, double B) {

    double cutoff = 1.0;

    if (q<1/B) return (    grad * grad / (1 - 1/B)     < cutoff);
    else       return (B * grad * grad / (q * (1 - q)) < cutoff);

}


// Evaluate whether or not a correlation is close to target value

bool isCloseC(double pij, double pi, double pj, double grad, double B) {

    double cutoff = 1.0;
    
    if (pij < 1/B) pij = 1/B;
    if ( pi < 1/B)  pi = 1/B;
    if ( pj < 1/B)  pj = 1/B;
    
    double dpij = sqrt(pij * (1 - pij) / B);
    double dpi  = sqrt( pi * (1 -  pi) / B);
    double dpj  = sqrt( pj * (1 -  pj) / B);
    
    double denom = dpij + pi * dpj + pj * dpi;

    return (grad * grad / (denom * denom) < cutoff);

}


// Chop off excessively large couplings/fields

void chop(Vector &J) {

    double cutoff = 15.0;

    for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) {
    
        if (J[i][j]<-1*cutoff) J[i][j] = -1 * cutoff;
        if (J[i][j]>cutoff)    J[i][j] = cutoff;
        
    } }

}


// Runs the main program

int runLearn(RunParametersQLS &r) {

    // Define variables
    
    Vector q, p;
    Vector J, expJ;
    
    // Get reference correlations from file
    
    FILE *compIn=fopen(r.getCompareInfile().c_str(),"r");
        
    if (compIn!=NULL) getCorrelations(compIn,q);
    else { printf("Error reading input from file %s",r.getCompareInfile().c_str()); return EXIT_FAILURE; }
    fclose(compIn);
    
    p.resize(q.size(),std::vector<double>());
    for (int i=0;i<p.size();i++) p[i].resize(q[i].size(), 0);
    
    // Retrieve initial couplings from file
    
    FILE *dataIn=fopen(r.getInfile().c_str(),"r");
        
    if (dataIn!=NULL) getCouplings(dataIn,J);
    else { printf("Error reading input from file %s",r.getInfile().c_str()); return EXIT_FAILURE; }
    fclose(dataIn);
    
    FILE *compOut=fopen(r.getCompareOutfile().c_str(),"w");
    
    // Resize expJ
    
    for (int i=0;i<J.size();i++) expJ.push_back(std::vector<double>(J[i].size(),0));
    for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) expJ[i][j] = exp(J[i][j]); }
    
    if (r.useVerbose) printf("Running sparse learning algorithm with target eps = %.2e...\n",r.epsilon);
    
    // Choose gauge
    
    if (r.useGI) {
    
        updateStep_ptr=&updateStep_GI;
        epsilonP2_ptr=&epsilonP2_GI;
        epsilonC_ptr=&epsilonC_GI;
        getMaxError_ptr=&getMaxError_GI;
        
    }
    else {
    
        updateStep_ptr=&updateStep;
        epsilonP2_ptr=&epsilonP2;
        epsilonC_ptr=&epsilonC;
        getMaxError_ptr=&getMaxError;
        
    }
	
    // Prepare Monte Carlo
    
    int    N        = sizetolength(J.size());   // System size
    double epsMax   = 1;                        // Maximum single term error tolerance (relative to maximum error norm, see below)
    double gamma    = 0;                        // Regularization strength (L2, set below)
    
    if (r.useGamma) {
    
        if (r.gamma==0) gamma=computeGamma_L2(q,r.sampleB);
        else            gamma=r.gamma;
        
    }
    
    if (r.useVerbose) printf("Got N=%d, len(h[0])=%d\n",N,(int)J[0].size());
    
    Vector grad(p);
    Vector weight(p);
    
    for (int i=0;i<p.size();i++) { for (int j=0;j<p[i].size();j++) { grad[i][j]=0; weight[i][j]=1; } }
    
    // Choose couplings to update
    
    IntVector nz;   // Only "nonzero" couplings will be updated when the learning algorithm is run
    
    for (int i=0;i<N;i++) { for (int j=i+1;j<N;j++) {
        
        // Add all pairs (future, update only select couplings)
        
        std::vector<int> pair(2,0);
        pair[0]=i;
        pair[1]=j;
            
        nz.push_back(pair);
            
    } }
    
    // Get default starting configuration, if nontrivial
    
    std::vector<int> lattice(N);
    
    if (r.useStart) {
    
        FILE *startIn=fopen(r.getStartInfile().c_str(),"r");
        for (int i=0;i<N;i++) fscanf(startIn,"%d",&lattice[i]);
        fclose(startIn);
    
    }
    
    else { for (int i=0;i<N;i++) lattice[i]=(int) q[i].size(); }
    
    // Prepare to simulate
    
    srand((unsigned)time(0));
    
    // Run MC and get error
    
    int count=0;
    std::vector<double> error(3,0);
    
    getErrorMCLearn(q, J, expJ, r.sampleB, r.b, r.runs, gamma, alpha, p, error, lattice);
    
    // Loop Monte Carlo until convergence
    
    if (r.useVerbose) printf("%d\t%.6e\t%.6e\t%.6e\t%.6e\n",count,error[0],error[1],error[2],LInfinity(weight));
    
    while ( error[0]>r.epsilon || error[1]>r.epsilon || error[2]>r.epsilon ) {
        
        // Make step and update gradient and weights
        
        (*updateStep_ptr)(q, p, nz, gamma, r.sampleB, J, expJ, grad, weight);
        count++;
        
        // Run MC and get error
        
        std::vector<double> lastError(error);
        
        getErrorMCLearn(q, J, expJ, r.sampleB, r.b, r.runs, gamma, alpha, p, error, lattice);
        
        // Output progress
        
        if (r.useVerbose) printf("%d\t%.6e\t%.6e\t%.6e\t%.6e\n",count,error[0],error[1],error[2],LInfinity(weight));
        fprintf(compOut,"%d\t%.6e\t%.6e\t%.6e\t%.6e\n",count,error[0],error[1],error[2],LInfinity(weight));
        fflush(compOut);
        
        // Write out couplings
        
        FILE *dataOut=fopen(r.getCouplingsOutfile().c_str(),"w");
        
        if (dataOut!=NULL) printCouplings(dataOut, J);
        else { printf("Error writing output to file %s",r.getCouplingsOutfile().c_str()); return EXIT_FAILURE;; }
        fclose(dataOut);
    
    }
    
    // Write out couplings
    
    FILE *dataOut=fopen(r.getCouplingsOutfile().c_str(),"w");
    
    if (dataOut!=NULL) printCouplings(dataOut, J);
    else { printf("Error writing output to file %s",r.getCouplingsOutfile().c_str()); return EXIT_FAILURE;; }
    fclose(dataOut);
    
    return EXIT_SUCCESS;
    
}

