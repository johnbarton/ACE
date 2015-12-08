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

#include "tools.h"      // Numerical tools
#include "monteCarlo.h" // MC declarations
#include "mtrnd.h"      // Random number generator


// GLOBAL VARIABLES

double (*epsilonP2_ptr)(const Vector &, const Vector &, int, double, const Vector &, double);
double (*epsilonC_ptr)(const Vector &, const Vector &, int, double, const Vector &, double);
double (*getMaxError_ptr)(const Vector &, const Vector &, double, const Vector &, double, double);




// Compute epsilonP for two vectors of correlations

double epsilonP(const Vector &q, const Vector &p, int N, double maxPrecision, const Vector &J, double gamma, double alpha) {
    
    double err  = 0.0;
    double Neff = 0.0;
    
    for (int i=0;i<N;i++) { for (int a=0;a<q[i].size();a++) {
    
        Neff++;
        double gradError = pow(q[i][a] - p[i][a] - (2 * gamma * alpha * J[i][a]), 2.0);
    
        if (q[i][a]<maxPrecision) err += gradError / (maxPrecision * (1 - maxPrecision));
        else                      err += gradError / (q[i][a]      * (1 -      q[i][a]));
        
    } }
    
    return sqrt( err / (Neff * maxPrecision) );
    
}


// Compute epsilonP2 for two vectors of correlations

double epsilonP2(const Vector &q, const Vector &p, int N, double maxPrecision, const Vector &J, double gamma) {
    
    double err   = 0;
    double NJeff = 0;
    
    for (int i=0;i<N;i++) { for (int j=i+1;j<N;j++) {
    
        int idx = index(i,j,N);
        
        for (int a=0;a<J[idx].size();a++) {
        
            NJeff++;

            double K         = J[idx][a];
            double gradError = pow(p[idx][a] - q[idx][a] + (2 * gamma * K), 2.0);
    
            if (q[idx][a]<maxPrecision) err += gradError/(maxPrecision * (1 - maxPrecision));
            else                        err += gradError/(   q[idx][a] * (1 -    q[idx][a]));

        } }
            
    }
    
    return sqrt( err / (NJeff * maxPrecision) );
    
}


// Compute epsilonP2 for two vectors of correlations

double epsilonP2_GI(const Vector &q, const Vector &p, int N, double maxPrecision, const Vector &J, double gamma) {
    
    double err   = 0;
    double NJeff = 0;
    
    for (int i=0;i<N;i++) { for (int j=i+1;j<N;j++) {
    
        int idx = index(i,j,N);
        
        //Get quantities for computing K_ij^ab
    
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
        
            NJeff++;

            int    sab       = sindex(a,b,J[i].size(),J[j].size());
            double K         = J[idx][sab] - (suma[a] * inv_qj) - (sumb[b] * inv_qi) + (sumall * inv_qi * inv_qj);
            double gradError = pow(p[idx][sab] - q[idx][sab] + (2 * gamma * K), 2.0);
    
            if (q[idx][sab]<maxPrecision) err += gradError/(maxPrecision * (1 - maxPrecision));
            else                          err += gradError/( q[idx][sab] * (1 - q[idx][sab]));

        } }
            
    } }
    
    return sqrt( err / (NJeff * maxPrecision) );
    
}


// Compute epsilonC for two vectors of correlations

double epsilonC(const Vector &q, const Vector &p, int N, double maxPrecision, const Vector &J, double gamma) {
    
    double err   = 0;
    double NJeff = 0;
    
    for (int i=0;i<N;i++) { for (int j=i+1;j<N;j++) {
    
        int idx = index(i,j,N);
        
        for (int a=0;a<J[i].size();a++) { for (int b=0;b<J[j].size();b++) {
        
            NJeff++;

            int    sab       = sindex(a,b,J[i].size(),J[j].size());
            double K         = J[idx][sab];
            double gradError = pow((q[idx][sab]-q[i][a]*q[j][b]) - (p[idx][sab]-p[i][a]*p[j][b]) - (2 * gamma * K) , 2.0);
            
            if (q[idx][sab]<maxPrecision) err += gradError / pow(sqrt(maxPrecision * (1 - maxPrecision)) + q[i][a] * sqrt(q[j][b] * (1 - q[j][b])) + q[j][b] * sqrt(q[i][a] * (1 - q[i][a])),2);
            else                          err += gradError / pow(sqrt(q[idx][sab]  * (1 -  q[idx][sab])) + q[i][a] * sqrt(q[j][b] * (1 - q[j][b])) + q[j][b] * sqrt(q[i][a] * (1 - q[i][a])),2);

        } }
            
    } }
    
    return sqrt( err / (NJeff * maxPrecision) );
    
}


// Compute epsilonC for two vectors of correlations

double epsilonC_GI(const Vector &q, const Vector &p, int N, double maxPrecision, const Vector &J, double gamma) {
    
    double err   = 0;
    double NJeff = 0;
    
    for (int i=0;i<N;i++) { for (int j=i+1;j<N;j++) {
    
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
        
            NJeff++;

            int    sab       = sindex(a,b,J[i].size(),J[j].size());
            double K         = J[idx][sab] - (suma[a] * inv_qj) - (sumb[b] * inv_qi) + (sumall * inv_qi * inv_qj);
            double gradError = pow((q[idx][sab]-q[i][a]*q[j][b]) - (p[idx][sab]-p[i][a]*p[j][b]) - (2 * gamma * K) , 2.0);
            
            if (q[idx][sab]<maxPrecision) err += gradError / pow(sqrt(maxPrecision * (1 - maxPrecision)) + q[i][a] * sqrt(q[j][b] * (1 - q[j][b])) + q[j][b] * sqrt(q[i][a] * (1 - q[i][a])),2);
            else                          err += gradError / pow(sqrt(q[idx][sab]  * (1 -  q[idx][sab])) + q[i][a] * sqrt(q[j][b] * (1 - q[j][b])) + q[j][b] * sqrt(q[i][a] * (1 - q[i][a])),2);

        } }
            
    } }
    
    return sqrt( err / (NJeff * maxPrecision) );
    
}


// Compute maximum error on the one- and two-point correlations (sparse)

double getMaxError(const Vector &q, const Vector &p, const IntVector &nz, double maxPrecision) {

    int    N        = sizetolength(q.size());
    double NJeff    = 0;
    double maxError = 0;
    
    // Check all single sites
        
    for (int i=0;i<N;i++) { for (int j=0;j<q[i].size();j++) {
    
        double gradError = pow(q[i][j] - p[i][j], 2.0);
        
        if (q[i][j]<maxPrecision) gradError /= (1 - maxPrecision);
        else                      gradError /= (maxPrecision * q[i][j] * (1 - q[i][j]));
            
        if (gradError>maxError) maxError=gradError;
        
    } }
        
    // Check nonzero pairs
        
    for (int n=0;n<nz.size();n++) {
    
        int i=nz[n][0];
        int j=nz[n][1];
            
        int idx=index(i,j,N);
            
        for (int a=0;a<q[idx].size();a++) {
        
            NJeff++;
            
            double gradError = pow(q[idx][a] - p[idx][a], 2.0);
            
            if (q[idx][a]<maxPrecision) gradError /= (1 - maxPrecision);
            else                        gradError /= (maxPrecision * q[idx][a] * (1 - q[idx][a]));

            if (gradError>maxError) maxError=gradError;
            
        }
            
    }
    
    return sqrt( maxError / (2 * log(NJeff) ) );

}


// Compute maximum error on the one- and two-point correlations (dense, including regularization)

double getMaxError(const Vector &q, const Vector &p, double maxPrecision, const Vector &J, double gamma, double alpha) {

    int    N        = sizetolength(q.size());
    double NJeff    = 0;
    double maxError = 0;
    
    // Check all single sites
        
    for (int i=0;i<N;i++) { for (int j=0;j<q[i].size();j++) {
    
        double gradError = pow(q[i][j] - p[i][j] - (2 * gamma * alpha * J[i][j]) , 2.0);
        
        if (q[i][j]<maxPrecision) gradError /= (1 - maxPrecision);
        else                      gradError /= (maxPrecision * q[i][j] * (1 - q[i][j]));
            
        if (gradError>maxError) maxError = gradError;
        
    } }
        
    // Check all pairs
        
    for (int i=0;i<N;i++) { for (int j=i+1;j<N;j++) {
    
        int idx = index(i,j,N);
        
        for (int a=0;a<J[idx].size();a++) {
        
            NJeff++;

            double K         = J[idx][a];
            double gradError = pow(p[idx][a] - q[idx][a] + (2 * gamma * K), 2.0);
    
            if (q[idx][a]<maxPrecision) gradError /= (1 - maxPrecision);
            else                        gradError /= (maxPrecision * q[idx][a] * (1 - q[idx][a]));

            if (gradError>maxError) maxError=gradError;
            
        }
            
    } }
    
    return sqrt( maxError / sqrt(2 * log(NJeff) ) );

}


// Compute maximum error on the one- and two-point correlations (dense, including regularization)

double getMaxError_GI(const Vector &q, const Vector &p, double maxPrecision, const Vector &J, double gamma, double alpha) {

    int    N        = sizetolength(q.size());
    double NJeff    = 0;
    double maxError = 0;
    
    // Check all single sites
        
    for (int i=0;i<N;i++) { for (int j=0;j<q[i].size();j++) {
    
        double gradError = pow(q[i][j] - p[i][j] - (2 * gamma * alpha * J[i][j]) , 2.0);
        
        if (q[i][j]<maxPrecision) gradError /= (1 - maxPrecision);
        else                      gradError /= (maxPrecision * q[i][j] * (1 - q[i][j]));
            
        if (gradError>maxError) maxError = gradError;
        
    } }
        
    // Check all pairs
        
    for (int i=0;i<N;i++) { for (int j=i+1;j<N;j++) {
    
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
        
            NJeff++;

            int    sab       = sindex(a,b,J[i].size(),J[j].size());
            double K         = J[idx][sab] - (suma[a] * inv_qj) - (sumb[b] * inv_qi) + (sumall * inv_qi * inv_qj);
            double gradError = pow(p[idx][sab] - q[idx][sab] + (2 * gamma * K), 2.0);
    
            if (q[idx][sab]<maxPrecision) gradError /= (1 - maxPrecision);
            else                          gradError /= (maxPrecision * q[idx][sab] * (1 - q[idx][sab]));

            if (gradError>maxError) maxError=gradError;
            
        } }
            
    } }
    
    return sqrt( maxError / sqrt(2 * log(NJeff) ) );

}


// Set energies based on an initial configuration

void setEnergies(const Vector &expJ, const IntVector &neighbor, const std::vector<int> &lattice, const std::vector<int> &nonzero, Vector &expH, std::vector<double> &sumExpH) {
    
    // Get energy of each state at each site (modulo fields/couplings at other sites)
    
    for (int i=0;i<sumExpH.size();i++) { sumExpH[i] = 1; for (int j=0;j<expH[i].size();j++) expH[i][j]=1; }

    for (int i=0;i<expH.size();i++) { for (int j=0;j<expH[i].size()-1;j++) {
        
        expH[i][j] = expJ[i][j];
            
        for (int k=1;k<neighbor[i].size();k++) { if (nonzero[neighbor[i][k]]) {
        
            if (i<neighbor[i][k]) expH[i][j] *= expJ[index(i,neighbor[i][k],lattice.size())][sindex(j,lattice[neighbor[i][k]],expJ[i].size(),expJ[neighbor[i][k]].size())];
            else                  expH[i][j] *= expJ[index(neighbor[i][k],i,lattice.size())][sindex(lattice[neighbor[i][k]],j,expJ[neighbor[i][k]].size(),expJ[i].size())];
                
        } }
            
        sumExpH[i] += expH[i][j];
            
    } }
    
}


// Update energies after a successful spin flip

void updateEnergies(const Vector &expJ, const IntVector &neighbor, const std::vector<int> &lattice, const std::vector<int> &nonzero, int x, int prev, Vector &expH, std::vector<double> &sumExpH) {

    // Initial state zero
    
    if (prev==expJ[x].size()) {
    
        for (int i=1;i<neighbor[x].size();i++) {
            
            int n = neighbor[x][i];
            
            sumExpH[n] = 1;
            
            if (x<n) { for (int j=0;j<expH[n].size()-1;j++) {
            
                expH[n][j] *= expJ[index(x,n,lattice.size())][sindex(lattice[x],j,expJ[x].size(),expJ[n].size())];
                sumExpH[n] += expH[n][j];
                
            } }
                
            else { for (int j=0;j<expH[n].size()-1;j++) {
            
                expH[n][j] *= expJ[index(n,x,lattice.size())][sindex(j,lattice[x],expJ[n].size(),expJ[x].size())];
                sumExpH[n] += expH[n][j];
                
            } }
            
        }
        
    }
    
    // Final state zero
    
    else if (lattice[x]==expJ[x].size()) {
    
        for (int i=1;i<neighbor[x].size();i++) {
            
            int n = neighbor[x][i];
            
            sumExpH[n] = 1;
            
            if (x<n) { for (int j=0;j<expH[n].size()-1;j++) {
                
                expH[n][j] /= expJ[index(x,n,lattice.size())][sindex(prev,j,expJ[x].size(),expJ[n].size())];
                sumExpH[n] += expH[n][j];
                    
            } }
            
            else { for (int j=0;j<expH[n].size()-1;j++) {
                
                expH[n][j] /= expJ[index(n,x,lattice.size())][sindex(j,prev,expJ[n].size(),expJ[x].size())];
                sumExpH[n] += expH[n][j];
                
            } }
        
        }
        
    }
    
    // Both initial and final state are nonzero
    
    else {
    
        for (int i=1;i<neighbor[x].size();i++) {
            
            int n = neighbor[x][i];
            
            sumExpH[n] = 1;
            
            if (x<n) { for (int j=0;j<expH[n].size()-1;j++) {
                
                expH[n][j] *= expJ[index(x,n,lattice.size())][sindex(lattice[x],j,expJ[x].size(),expJ[n].size())]/expJ[index(x,n,lattice.size())][sindex(prev,j,expJ[x].size(),expJ[n].size())];
                sumExpH[n] += expH[n][j];
                
            } }
            
            else { for (int j=0;j<expH[n].size()-1;j++) {
            
                expH[n][j] *= expJ[index(n,x,lattice.size())][sindex(j,lattice[x],expJ[n].size(),expJ[x].size())]/expJ[index(n,x,lattice.size())][sindex(j,prev,expJ[n].size(),expJ[x].size())];
                sumExpH[n] += expH[n][j];
                
            } }
            
        }
        
    }

}


// Rejection-free Monte Carlo step (update energies on flip only)

void dynamicsRF(const Vector &expJ, const IntVector &neighbor, int x, double y, std::vector<int> &lattice, std::vector<int> &nonzero, Vector &expH, std::vector<double> &sumExpH) {
    
    // Check for spin flip, choose new state
    
    if ( y < ((sumExpH[x] - expH[x][lattice[x]]) / sumExpH[x]) ) {
        
        double tempSum=0;
        int prev=lattice[x];
        
        for (int i=0;i<expH[x].size();i++) { if (i!=lattice[x]) {
            
            tempSum += expH[x][i]/sumExpH[x];
            
            if (y<tempSum) {
                
                if (lattice[x]==expJ[x].size()) nonzero[x]=1;
                else if (i==expJ[x].size())     nonzero[x]=0;
                
                lattice[x]=i;
                break;
                
            }
            
        } }
        
        updateEnergies(expJ, neighbor, lattice, nonzero, x, prev, expH, sumExpH);
        
    }

}


// Update correlations based on the current lattice configuration

void updateCorrelations(const std::vector<int> &lattice, const std::vector<int> &nonzero, Vector &p) {

    for (int i=0;i<lattice.size();i++) { if (nonzero[i]) {
            
        p[i][lattice[i]]+=1;
            
        int off=offset(i, lattice.size());
            
        for (int j=i+1;j<lattice.size();j++) { if (nonzero[j]) {
        
            p[off + j][sindex(lattice[i],lattice[j],p[i].size(),p[j].size())]+=1;
                
        } }
            
    } }

}


// Update correlations for Generative Tests based on the current lattice configuration with 3-point correlations

void updateCorrelations(const std::vector<int> &lattice, const std::vector<int> &nonzero, Vector &p, std::vector<double> &pk, std::vector<std::vector<std::vector<std::vector<double> > > > &p3, std::vector<int> &cons) {
    
    for (int i=0;i<lattice.size();i++) { if (nonzero[i]) {
        
        p[i][lattice[i]]+=1;
            
        int off=offset(i, lattice.size());
            
        for (int j=i+1;j<lattice.size();j++) { if (nonzero[j]) {
        
            p[off + j][sindex(lattice[i],lattice[j],p[i].size(),p[j].size())]+=1;
	    
            for (int k=j+1;k<lattice.size();k++) {  if (nonzero[k]) {
				    
                p3[i][j][k][sindex3(lattice[i],lattice[j],lattice[k],p[i].size(),p[j].size(),p[k].size())]+=1;

            } }
                
        } }
        
    } }
    
    int nmut=0;
    for (int i=0;i<lattice.size();i++) {
    
        if (lattice[i]!=cons[i]) nmut++;
        else if (cons[i]>=0)     nmut++;
        
    }
    
    pk[nmut]+=1;

}


// Update correlations for Generative Tests based on the current lattice configuration WITHOUT 3-point correlations

void updateCorrelations(const std::vector<int> &lattice, const std::vector<int> &nonzero, Vector &p, std::vector<double> &pk, std::vector<int> &cons) {
    
    for (int i=0;i<lattice.size();i++) { if (nonzero[i]) {
        
        p[i][lattice[i]]+=1;
            
        int off=offset(i, lattice.size());
            
        for (int j=i+1;j<lattice.size();j++) { if (nonzero[j]) {
        
            p[off + j][sindex(lattice[i],lattice[j],p[i].size(),p[j].size())]+=1;
                
        } }
        
    } }
    
    int nmut=0;
    for (int i=0;i<lattice.size();i++) {
    
        if (lattice[i]!=cons[i]) nmut++;
        else if (cons[i]>=0)     nmut++;
        
    }
    
    pk[nmut]+=1;

}


// Get Monte Carlo sample

void getSample(const Vector &expJ, const IntVector &neighbor, MT::MersenneTwist &mt, int N, int numSteps, int stepSize, std::vector<int> &lattice, std::vector<int> &nonzero, Vector &p, Vector &expH, std::vector<double> &sumExpH) {

    for (int n=0;n<numSteps;n++) {
    
        //for (int k=0;k<stepSize;k++) dynamicsRF(expJ, neighbor, mt.genrand_int31()%N, mt.genrand_real1(), lattice, nonzero, expH, sumExpH);
        for (int k=0;k<stepSize;k++) dynamicsRF(expJ, neighbor, k%N, mt.genrand_real1(), lattice, nonzero, expH, sumExpH);
        updateCorrelations(lattice, nonzero, p);
        
    }

}


// Get Monte Carlo sample for GenTest with 3-point correlations

void getSampleGenTest(const Vector &expJ, const IntVector &neighbor, MT::MersenneTwist &mt, int N, int numSteps, int stepSize, std::vector<int> &lattice, std::vector<int> &nonzero, Vector &p, Vector &expH, std::vector<double> &sumExpH,  std::vector<double> &pk, std::vector<std::vector<std::vector<std::vector<double> > > > &p3,std::vector<int> &cons) {

    for (int n=0;n<numSteps;n++) {
    
        for (int k=0;k<stepSize;k++) dynamicsRF(expJ, neighbor, mt.genrand_int31()%N, mt.genrand_real1(), lattice, nonzero, expH, sumExpH);
        updateCorrelations(lattice, nonzero, p, pk, p3, cons);
    }

}

// Get Monte Carlo sample for GenTest WITHOUT 3-point correlations

void getSampleGenTest(const Vector &expJ, const IntVector &neighbor, MT::MersenneTwist &mt, int N, int numSteps, int stepSize, std::vector<int> &lattice, std::vector<int> &nonzero, Vector &p, Vector &expH, std::vector<double> &sumExpH,  std::vector<double> &pk,std::vector<int> &cons) {

    for (int n=0;n<numSteps;n++) {
    
        for (int k=0;k<stepSize;k++) dynamicsRF(expJ, neighbor, mt.genrand_int31()%N, mt.genrand_real1(), lattice, nonzero, expH, sumExpH);
        updateCorrelations(lattice, nonzero, p, pk, cons);
    }

}



// Get effective number of fields and couplings

void getNeff(const Vector &q, int N, int &Neff, int &NJeff) {

    for (int i=0;i<N;i++) {
    
        Neff+=(int) q[i].size();
        for (int j=i+1;j<N;j++) NJeff+=(int) (q[i].size()*q[j].size());
        
    }

}


// Get the set of interacting "neighbors" for each site

void getNeighbors(const Vector &J, int N, double threshold, IntVector &neighbor) {

    for (int i=0;i<N;i++) neighbor[i].push_back(i);
    
    for (int i=0;i<N;i++) {
        
        int off=offset(i, N);
        
        for (int j=i+1;j<N;j++) { for (int k=0;k<J[off + j].size();k++) {
            
            if (fabs(J[off + j][k])>threshold) {
                
                neighbor[i].push_back(j);
                neighbor[j].push_back(i);
                break;
               
            }
                
        } }
        
    }

}


// Compute the inferrence error for a given set of couplings using rejection-free Monte Carlo (getError_RF)

void getError(const Vector &q, const Vector &J, int N, double B, int numSteps, int numRuns, double gamma, double alpha, std::vector<double> &error) {

    // Create and initialize Monte Carlo variables
    
    std::vector<double> lastError(2,0);
    Vector p(q.size(),std::vector<double>());
    for (int i=0;i<p.size();i++) p[i].resize(q[i].size(),0);
    
    std::vector<int> lattice(N,0);
    std::vector<int> nonzero(N,0);
    
    // Couplings and neighbors
    
    IntVector neighbor(N);
    getNeighbors(J, N, 0, neighbor);
    Vector expJ(J);
    for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) expJ[i][j]=exp(J[i][j]); }
    
    // Energies
    
    Vector expH;
    std::vector<double> sumExpH;
    
    sumExpH.resize(lattice.size(),0);
    expH.resize(lattice.size(),std::vector<double>());
    for (int i=0;i<expH.size();i++) expH[i].resize(expJ[i].size()+1,1);
    
    // MC variables
    
    double multiplier   = 1;
    double multMax      = 64;
    int    stepSize     = 4*N;
    double maxPrecision = 1/B;
    bool   runMC        = true;
    
    // Set variables for parallel
    
    int numThreads=1;
    
    srand((unsigned)time(0));
    unsigned long seed[numThreads];
    for (int i=0;i<numThreads;i++) seed[i]=rand();
    
    for (int i=0;i<lattice.size();i++) { lattice[i]=(int) q[i].size(); nonzero[i]=0; }

    MT::MersenneTwist mt;
    mt.init_genrand(seed[0]);
    
    // Thermalize and get initial error
    
    setEnergies(expJ, neighbor, lattice, nonzero, expH, sumExpH);
    getSample(expJ, neighbor, mt, N, numSteps, stepSize, lattice, nonzero, p, expH, sumExpH);
    
    // Normalize correlations and compute initial error

    for (int i=0;i<q.size();i++) { for (int j=0;j<q[i].size();j++) p[i][j] /= (double) (numThreads * numSteps * multiplier); }
    
    error[0] = epsilonP(       q, p, N, maxPrecision, J, gamma, alpha);
    error[1] = (*epsilonC_ptr)(q, p, N, maxPrecision, J, gamma);
    error[2] = (*getMaxError_ptr)(q, p, maxPrecision, J, gamma, alpha);
    
    for (int i=0;i<q.size();i++) { for (int j=0;j<q[i].size();j++) p[i][j] = 0; }
    
    // Go to MC
    
    while (runMC) {
                            
        for (int i=0;i<lattice.size();i++) { lattice[i]=(int) q[i].size(); nonzero[i] = 0; }
    
        setEnergies(expJ, neighbor, lattice, nonzero, expH, sumExpH);
        getSample(expJ, neighbor, mt, N, numSteps*multiplier, stepSize, lattice, nonzero, p, expH, sumExpH);
    
        // Normalize correlations and compute error

        for (int i=0;i<q.size();i++) { for (int j=0;j<q[i].size();j++) p[i][j] /= (double) (numThreads * numSteps * multiplier); }
    
        lastError[0] = error[0];
        lastError[1] = error[1];
    
        error[0] = epsilonP(       q, p, N, maxPrecision, J, gamma, alpha);
        error[1] = (*epsilonC_ptr)(q, p, N, maxPrecision, J, gamma);
        error[2] = (*getMaxError_ptr)(q, p, maxPrecision, J, gamma, alpha);
    
        // Check exit condition
            
        double deltaP = fabs(error[0] - lastError[0]);
        double deltaC = fabs(error[1] - lastError[1]);
  
        bool pBad  = ( (error[0] > 1.0) && (error[0] - 3 * deltaP > 1) );
        bool cBad  = ( (error[1] > 1.0) && (error[1] - 3 * deltaC > 1) );
        bool pGood = ( (error[0] < 1.0) && (error[0] + 1 * deltaP < 1) );
        bool cGood = ( (error[1] < 1.0) && (error[1] + 1 * deltaC < 1) );
            
        if ( pBad || cBad || (pGood && cGood) || multiplier==multMax ) runMC=false;
        
        else {
            
            for (int i=0;i<q.size();i++) { for (int j=0;j<q[i].size();j++) p[i][j] = 0; }
            multiplier*=2;
            //printf("m=%d, e_p=%.4e, e_c=%.4e\n",(int)multiplier,error[0],error[1]);
                
        }
    
    }

}


// Compute the inferrence error for a given set of couplings using rejection-free Monte Carlo (getError_RF)

void getErrorMCLearn(const Vector &q, const Vector &J, Vector &expJ, double B, int numSteps, int numRuns, double gamma, double alpha, Vector &p, std::vector<double> &error, std::vector<int> &latticeStart) {

    // Create and initialize Monte Carlo variables

    for (int i=0;i<p.size();i++) { for (int j=0;j<p[i].size();j++) p[i][j]=0; }

    int N = sizetolength(J.size());
    
    std::vector<int> lattice(N,0);
    std::vector<int> nonzero(N,0);
    
    // Neighbors
    
    IntVector neighbor(N);
    getNeighbors(J, N, 1.0e-5, neighbor);
    
    // Energies
    
    Vector expH;
    std::vector<double> sumExpH;
    
    sumExpH.resize(lattice.size(),0);
    expH.resize(lattice.size(),std::vector<double>());
    for (int i=0;i<expH.size();i++) expH[i].resize(expJ[i].size()+1,1);
    
    // MC variables

    int    stepSize     = 4*N;
    int    mutCut       = (N/2 < 80) ? N/2 : 120;
    double maxPrecision = 1/B;
    
    // Set variables for parallel
    
    int numThreads=1;
    
    srand((unsigned)time(0));
    unsigned long seed[numThreads];
    for (int i=0;i<numThreads;i++) seed[i]=rand();

    MT::MersenneTwist mt;
    mt.init_genrand(seed[0]);
        
    // Thermalize and get initial error
        
    for (int m=0;m<numRuns;m++) {
        
        for (int i=0;i<lattice.size();i++) {
        
            lattice[i]=latticeStart[i];
            if (lattice[i]==q[i].size()) nonzero[i]=0;
            else                         nonzero[i]=1;
            
        }
        
        setEnergies(expJ, neighbor, lattice, nonzero, expH, sumExpH);
        getSample(expJ, neighbor, mt, N, numSteps, stepSize, lattice, nonzero, p, expH, sumExpH);
            
    }
        
    // Normalize correlations and compute error
    
    for (int i=0;i<q.size();i++) { for (int j=0;j<q[i].size();j++) p[i][j] /= (double) (numThreads * numSteps * numRuns); }
        
    error[0] = epsilonP(        q, p, N, maxPrecision, J, gamma, alpha);
    error[1] = (*epsilonP2_ptr)(q, p, N, maxPrecision, J, gamma);
    error[2] = (*getMaxError_ptr)( q, p, maxPrecision, J, gamma, alpha);

}


// Compute the inferrence error for a given set of couplings using rejection-free Monte Carlo (getError_RF) with 3-point correlations

void getErrorGenTest(const Vector &J, Vector &expJ, double B, int numSteps, int numRuns, Vector &p, std::vector<int> &latticeStart, std::vector<double> &pk, std::vector<std::vector<std::vector<std::vector<double> > > > &p3, std::vector<int> &cons) {
    
    // Create and initialize Monte Carlo variables

    for (int i=0;i<p.size();i++) { for (int j=0;j<p[i].size();j++) p[i][j]=0; }

    int N = sizetolength(J.size());
    
    std::vector<int> lattice(N,0);
    std::vector<int> nonzero(N,0);
    
    // Neighbors
    
    IntVector neighbor(N);
    getNeighbors(J, N, 1.0e-5, neighbor);
    
    // Energies
    
    Vector expH;
    std::vector<double> sumExpH;
    
    sumExpH.resize(lattice.size(),0);
    expH.resize(lattice.size(),std::vector<double>());
    for (int i=0;i<expH.size();i++) expH[i].resize(expJ[i].size()+1,1);
    
    // MC variables

    int    stepSize     = 4*N;
    double maxPrecision = 1/B;
    
    // Set variables for parallel
    
    int numThreads=1;
    
    srand((unsigned)time(0));
    unsigned long seed[numThreads];
    for (int i=0;i<numThreads;i++) seed[i]=rand();

    MT::MersenneTwist mt;
    mt.init_genrand(seed[0]);
        
    // Thermalize and get initial error
        
    for (int m=0;m<numRuns;m++) {
        
        for (int i=0;i<lattice.size();i++) {
        
            lattice[i]=latticeStart[i];
            if (lattice[i]==p[i].size()) nonzero[i]=0;
            else                         nonzero[i]=1;
            
        }
        
        setEnergies(expJ, neighbor, lattice, nonzero, expH, sumExpH);
        getSampleGenTest(expJ, neighbor, mt, N, numSteps, stepSize, lattice, nonzero, p, expH, sumExpH, pk, p3, cons);
            
    }
        
    // Normalize correlations and compute error
    
    double Neff  = 0;
    double NJeff = 0;
    
    for (int i=0;i<lattice.size();i++) { for (int a=0;a<p[i].size();a++) {
	    
        Neff++;
        
        p[i][a]   /= (double) (numThreads * numSteps * numRuns);
	    
        for (int j=i+1;j<lattice.size();j++) { for (int b=0;b<p[j].size();b++) {
        
            NJeff++;
        
            int idx = index(i,j,N);
            int sab = sindex(a,b,p[i].size(),p[j].size());
		    
            p[idx][sab] /= (double) (numThreads * numSteps * numRuns);
		    
            for (int k=j+1;k<lattice.size();k++) { for (int c=0;c<p[k].size();c++) {
			    
                p3[i][j][k][sindex3(a,b,c,p[i].size(),p[j].size(),p[k].size())] /= (double) (numThreads * numSteps * numRuns);
                    
            } }
		
        } }
        
	} }

    for (int i=0;i<pk.size();i++) pk[i] /= (double) (numThreads * numSteps * numRuns);
    
}


// Compute the inferrence error for a given set of couplings using rejection-free Monte Carlo (getError_RF) WITHOUT 3-point correlations

void getErrorGenTest(const Vector &J, Vector &expJ, double B, int numSteps, int numRuns, Vector &p, std::vector<int> &latticeStart, std::vector<double> &pk, std::vector<int> &cons) {
    
    // Create and initialize Monte Carlo variables

    for (int i=0;i<p.size();i++) { for (int j=0;j<p[i].size();j++) p[i][j]=0; }

    int N = sizetolength(J.size());
    
    std::vector<int> lattice(N,0);
    std::vector<int> nonzero(N,0);
    
    // Neighbors
    
    IntVector neighbor(N);
    getNeighbors(J, N, 1.0e-5, neighbor);
    
    // Energies
    
    Vector expH;
    std::vector<double> sumExpH;
    
    sumExpH.resize(lattice.size(),0);
    expH.resize(lattice.size(),std::vector<double>());
    for (int i=0;i<expH.size();i++) expH[i].resize(expJ[i].size()+1,1);
    
    // MC variables

    int    stepSize     = 4*N;
    double maxPrecision = 1/B;
    
    // Set variables for parallel
    
    int numThreads=1;
    
    srand((unsigned)time(0));
    unsigned long seed[numThreads];
    for (int i=0;i<numThreads;i++) seed[i]=rand();

    MT::MersenneTwist mt;
    mt.init_genrand(seed[0]);
        
    // Thermalize and get initial error
        
    for (int m=0;m<numRuns;m++) {
        
        for (int i=0;i<lattice.size();i++) {
        
            lattice[i]=latticeStart[i];
            if (lattice[i]==p[i].size()) nonzero[i]=0;
            else                         nonzero[i]=1;
            
        }
        
        setEnergies(expJ, neighbor, lattice, nonzero, expH, sumExpH);
        getSampleGenTest(expJ, neighbor, mt, N, numSteps, stepSize, lattice, nonzero, p, expH, sumExpH, pk, cons);
            
    }
        
    // Normalize correlations and compute error
    
    double Neff  = 0;
    double NJeff = 0;
    
    for (int i=0;i<lattice.size();i++) { for (int a=0;a<p[i].size();a++) {
	    
        Neff++;
        
        p[i][a]   /= (double) (numThreads * numSteps * numRuns);
	    
        for (int j=i+1;j<lattice.size();j++) { for (int b=0;b<p[j].size();b++) {
        
            NJeff++;
        
            int idx = index(i,j,N);
            int sab = sindex(a,b,p[i].size(),p[j].size());
		    
            p[idx][sab] /= (double) (numThreads * numSteps * numRuns);
		
        } }
        
	} }

    for (int i=0;i<pk.size();i++) pk[i] /= (double) (numThreads * numSteps * numRuns);
    
}

