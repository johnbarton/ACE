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
#include <time.h>

#include "qGenTest.h"       // Generative Test declarations
#include "monteCarlo.h"     // Monte Carlo
#include "tools.h"          // Numerical tools
#include "io.h"             // Data I/O
#include "mtrnd.h"          // Random number generator





// Runs the main program

int runGenTest(RunParameters &r) {

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

    if (FILE *consIn = fopen(r.getConsensusInfile().c_str(),"r")) {
        getConsensus(consIn,cons);
        fclose(consIn);
    }
    else {
        printf("Problem retrieving data from file %s! The file may not exist, or it may be inaccessible\n",r.getConsensusInfile().c_str());
        return EXIT_FAILURE;
    }
    
    if (r.useVerbose) {

        printf("Reference sequence: ");
        for (int i=0;i<cons.size();i++) printf(" %d",cons[i]);
        printf("\n");

    }

    // Retrieve couplings from file
    
    if (FILE *dataIn=fopen(r.getInfile().c_str(),"r")) {
        getCouplings(dataIn,J);
        fclose(dataIn);
    }
    else {
        printf("Problem retrieving data from file %s! The file may not exist, or it may be inaccessible\n",r.getInfile().c_str());
        return EXIT_FAILURE;
    }
    

    // Resize expJ
    
    for (int i=0;i<J.size();i++) expJ.push_back(std::vector<double>(J[i].size(),0));
    for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) expJ[i][j] = exp(J[i][j]); }

    // Declare 2-point correlations, 3-point correlations, P(k) and magnetisations

    bool computeThreePoints = r.computeP3;
    int N = sizetolength(J.size()); // System size
    double gamma = 0;               // Regularization strength (L2, set below)
    
    if (r.useGamma) {

        if (r.gamma==0) gamma=1/(r.sampleB);
        else            gamma=r.gamma;
        
    }
    
    Vector p(J.size(),std::vector<double>());   // MC magnetisations and 2-point correlations
    Vector q(J.size(),std::vector<double>());   // MSA magnetisations and 2-point correlations
    Vector p3;                                  // MC 3-point correlations
    Vector q3;                                  // MSA 3-point correlations
    std::vector<double> pk(N+1,0);              // MC mutation probability
    std::vector<double> qk(N+1,0);              // MSA mutation probability
    
    for (int i=0;i<J.size();i++) {

        p[i].resize(J[i].size(),0);
        q[i].resize(J[i].size(),0);

    }
    
    if (computeThreePoints) {
        
        p3.resize(N*(N-1)*(N-2)/6, std::vector<double>());
        q3.resize(N*(N-1)*(N-2)/6, std::vector<double>());
        
        for (int i=0;i<N;i++) { for (int j=i+1;j<N;j++) { for (int k=j+1;k<N;k++) {

            p3[index(i,j,k,N)].resize(J[i].size()*J[j].size()*J[k].size(),0);
            q3[index(i,j,k,N)].resize(J[i].size()*J[j].size()*J[k].size(),0);
        
        } } }
        
    }

    // Get sequences from MSA file and compute correlations

    ;

    if (FILE *alIn = fopen(r.getInfileAl().c_str(),"r")) {
        FILE *weightIn = fopen(r.getWeights().c_str(),"r");
        
        if (computeThreePoints) getAlignment(alIn, weightIn, J, q, q3, qk, cons, r.MSAEnOutfile());
        else                    getAlignment(alIn, weightIn, J, q,     qk, cons, r.MSAEnOutfile());
        
        if (weightIn!=NULL) fclose(weightIn);
        fclose(alIn);
    }
    else {
        printf("Problem retrieving data from file %s! The file may not exist, or it may be inaccessible\n",r.getInfileAl().c_str());
        return EXIT_FAILURE;
    }

    if (r.useVerbose) printf("Got N=%d, len(h[0])=%d\n",N,(int)J[0].size());

    // Get default starting configuration, if nontrivial

    std::vector<int> lattice(N);

    if (r.useStart) {

        FILE *startIn=fopen(r.getStartInfile().c_str(),"r");
        for (int i=0;i<N;i++) fscanf(startIn,"%d",&lattice[i]);

    }

    else { for (int i=0;i<N;i++) lattice[i]=(int) p[i].size(); }

    // Run MC and get correlations
    
//    //DEBUG
//    for (int i=0;i<lattice.size();i++) printf("%d ",lattice[i]);
//    printf("\n");
//    for (int i=0;i<qk.size();i++) printf("%le ",qk[i]);
//    printf("\n");
//    for (int i=0;i<pk.size();i++) printf("%le ",pk[i]);
//    printf("\n");
    for (int i=0;i<cons.size();i++) { if (cons[i]<0) cons[i]=p[i].size(); }
//    //DEBUG

    srand((unsigned)time(0));

    if (computeThreePoints) getErrorGenTest(J, expJ, r.sampleB, r.b, r.runs, p, lattice, pk, p3, cons, r.recMSA, r.MSAOutfile(), r.EnergiesOutfile());
    else                    getErrorGenTest(J, expJ, r.sampleB, r.b, r.runs, p, lattice, pk,     cons, r.recMSA, r.MSAOutfile(), r.EnergiesOutfile());
    
    // Print epsilon errors before expanding
    
    double alpha        = 0.01;   // Regularization normalization for fields
    double maxPrecision = 1/(r.sampleB);
    double erelp1       = epsilonP(q, p, N, maxPrecision, J, gamma, alpha);
    double erelp2       = (*epsilonP2_ptr)(q, p, N, maxPrecision, J, gamma);
    double emax         = (*getMaxError_ptr)( q, p, maxPrecision, J, gamma, alpha);

    printf("Relative errors (w/out gauge correlations): P %f, P2 %f MAX %f gamma %f\n",erelp1,erelp2,emax,gamma);
    
    // Expand correlations to include gauge states
    
    if (computeThreePoints) {
        expandCorrelations3(q, q3);
        expandCorrelations3(p, p3);
    }
    else {
        Vector qgauge, pgauge;
        expandCorrelations(q, qgauge);
        expandCorrelations(p, pgauge);
        q = qgauge;
        p = pgauge;
    }

    // Compute and print correlations
    
    if (FILE *pkout = fopen(r.getPKOutfile().c_str(),"w")) {
        for (int i=0;i<qk.size();i++) fprintf(pkout,"%d %le %le\n",i,qk[i],pk[i]);
        fclose(pkout);
    }
    else {
        printf("Error writing output to file %s!\n",r.getPKOutfile().c_str());
        return EXIT_FAILURE;
    }
    
    FILE *p1out = fopen(r.getMOutfile().c_str(),"w");
    FILE *p2out = fopen(r.getP2Outfile().c_str(),"w");
    FILE *c2out = fopen(r.getCCOutfile().c_str(),"w");
    FILE *p3out = NULL;
    FILE *c3out = NULL;
    
    if (computeThreePoints) {
        p3out = fopen(r.getP3Outfile().c_str(),"w");
        c3out = fopen(r.getC3Outfile().c_str(),"w");
    }
    
    std::vector<double> sizes(5,0);             // denominator for averages
    for (int i=0;i<N;i++) {
        sizes[0] += p[i].size();
        for (int j=i+1;j<N;j++) {
            sizes[1] += p[i].size()*p[j].size();
            sizes[2] += p[i].size()*p[j].size();
            for (int k=j+1;k<N;k++) {
                sizes[3] += p[i].size()*p[j].size()*p[k].size();
                sizes[4] += p[i].size()*p[j].size()*p[k].size();
            }
        }
    }
    
    std::vector<std::vector<double> > averages; // hold correlation averages for computing R (p1, p2, c2, p3, c3)
    for (int i=0;i<5;i++) averages.push_back(std::vector<double>(5,0));
    
    double nmax       = (r.useNMax) ? r.nmax : sizes[3];    // rough maximum number of 3-point correlations to print
    double randthresh = nmax/sizes[3];                      // randomly select correlations for printing if random < randthresh OR size > pthresh
    double pthresh    = r.pthresh;                          // threshold correlation value for printing
    
    if (r.computePThresh) {     // default value: 10 x (avg q)^3
    
        double meanq = ((double) N) / sizes[0];
        pthresh      = 10.0 * pow(meanq, 3.0);
        
    }
    
    MT::MersenneTwist mt;
    mt.init_genrand(rand());
    
    // Loop through correlations
    
    for (int i=0;i<N;i++) { for (int a=0;a<p[i].size();a++) {

        averages[0][0] += q[i][a];              // average q1
        averages[0][1] += p[i][a];              // average p1
        averages[0][2] += q[i][a] * p[i][a];    // average q1 * p1
        averages[0][3] += q[i][a] * q[i][a];    // average q1 * q1
        averages[0][4] += p[i][a] * p[i][a];    // average p1 * p1
        
        fprintf(p1out,"%le\t%le\t%le\n",q[i][a],p[i][a],sqrt((maxPrecision*q[i][a]*(1-q[i][a]))));
        fflush(p1out);

        for (int j=i+1;j<N;j++) { for (int b=0;b<p[j].size();b++) {

            int idx = index(i,j,N);
            int sab = sindex(a,b,p[i].size(),p[j].size());
            
            averages[1][0] += q[idx][sab];                  // average q2
            averages[1][1] += p[idx][sab];                  // average p2
            averages[1][2] += q[idx][sab] * p[idx][sab];    // average q2 * p2
            averages[1][3] += q[idx][sab] * q[idx][sab];    // average q2 * q2
            averages[1][4] += p[idx][sab] * p[idx][sab];    // average p2 * p2
            
            fprintf(p2out,"%le\t%le\t%le\n",q[idx][sab],p[idx][sab],sqrt((maxPrecision*q[idx][sab]*(1-q[idx][sab]))));
            fflush(p2out);
            
            double qc2 = q[idx][sab] - (q[i][a] * q[j][b]);
            double pc2 = p[idx][sab] - (p[i][a] * p[j][b]);
            double ec2 = sqrt((maxPrecision*q[idx][sab]*(1-q[idx][sab])) + (q[i][a]*q[i][a]*(maxPrecision*q[j][b]*(1-q[j][b]))) + (q[j][b]*q[j][b]*(maxPrecision*q[i][a]*(1-q[i][a]))));
            
            averages[2][0] += qc2;          // average qc2
            averages[2][1] += pc2;          // average pc2
            averages[2][2] += qc2 * pc2;    // average qc2 * pc2
            averages[2][3] += qc2 * qc2;    // average qc2 * qc2
            averages[2][4] += pc2 * pc2;    // average pc2 * pc2
            
            fprintf(c2out,"%le\t%le\t%le\n",qc2,pc2,ec2);
            fflush(c2out);

            if (computeThreePoints) { for (int k=j+1;k<N;k++) { for (int c=0;c<p[k].size();c++) {

                int ijx  = idx;
                int ikx  = index(i,k,N);
                int jkx  = index(j,k,N);
                int sac  = sindex(a,c,p[i].size(),p[k].size());
                int sbc  = sindex(b,c,p[j].size(),p[k].size());
                int sabc = sindex3(a,b,c,p[i].size(),p[j].size(),p[k].size());
                
                double qp3 = q3[index(i,j,k,N)][sabc];
                double pp3 = p3[index(i,j,k,N)][sabc];
                
                averages[3][0] += qp3;          // average q3
                averages[3][1] += pp3;          // average p3
                averages[3][2] += qp3 * pp3;    // average q3 * p3
                averages[3][3] += qp3 * qp3;    // average q3 * q3
                averages[3][4] += pp3 * pp3;    // average p3 * p3
                
                if (qp3>pthresh || pp3>pthresh || mt.genrand_real1()<randthresh) {
                    fprintf(p3out,"%le\t%le\t%le\n",qp3,pp3,sqrt((maxPrecision*qp3*(1-qp3))));
                    fflush(p3out);
                }
                
                double qc3 = qp3 - (q[i][a]*q[jkx][sbc]) - (q[j][b]*q[ikx][sac]) - (q[k][c]*q[ijx][sab]) + (2*(q[i][a]*q[j][b]*q[k][c]));
                double pc3 = pp3 - (p[i][a]*p[jkx][sbc]) - (p[j][b]*p[ikx][sac]) - (p[k][c]*p[ijx][sab]) + (2*(p[i][a]*p[j][b]*p[k][c]));
                double ec3 = sqrt((maxPrecision*qp3*(1-qp3)) + (q[i][a]*q[i][a]*(maxPrecision*q[jkx][sbc]*(1-q[jkx][sbc]))) + (q[jkx][sbc]*q[jkx][sbc]*(maxPrecision*q[i][a]*(1-q[i][a]))) + (q[j][b]*q[j][b]*(maxPrecision*q[ikx][sac]*(1-q[ikx][sac]))) + (q[ikx][sac]*q[ikx][sac]*(maxPrecision*q[j][b]*(1-q[j][b]))) + (q[k][c]*q[k][c]*(maxPrecision*q[ijx][sab]*(1-q[ijx][sab]))) + (q[ijx][sab]*q[ijx][sab]*(maxPrecision*q[k][c]*(1-q[k][c]))) + (4*(((maxPrecision*q[i][a]*(1-q[i][a]))*q[j][b]*q[k][c]*q[j][b]*q[k][c]) + (q[i][a]*q[i][a]*(maxPrecision*q[j][b]*(1-q[j][b]))*q[k][c]*q[k][c]) + (q[i][a]*q[j][b]*q[i][a]*q[j][b]*(maxPrecision*q[k][c]*(1-q[k][c]))))));
            
                averages[4][0] += qc3;          // average qc3
                averages[4][1] += pc3;          // average pc3
                averages[4][2] += qc3 * pc3;    // average qc3 * pc3
                averages[4][3] += qc3 * qc3;    // average qc3 * qc3
                averages[4][4] += pc3 * pc3;    // average pc3 * pc3

                if (fabs(qc3)>pthresh || fabs(pc3)>pthresh || mt.genrand_real1()<randthresh) {
                    fprintf(c3out,"%le\t%le\t%le\n",qc3,pc3,ec3);
                    fflush(c3out);
                }

            } } }

        } }

    } }
    
    // Normalize averages
    
    for (int i=0;i<averages.size();i++) { for (int j=0;j<averages[i].size();j++) averages[i][j] /= sizes[i]; }
    
    // Compute and print errors
    
    std::vector<double> RMSErr(5,0);    // sqrt( q^2 + p^2 - 2 qp )
    for (int i=0;i<RMSErr.size();i++) RMSErr[i] = sqrt(averages[i][3] + averages[i][4] - (2*averages[i][2]));

    if (computeThreePoints) printf("Root mean square deviations: P %f, P2 %f, P3 %f, C2 %f, C3 %f \n",RMSErr[0],RMSErr[1],RMSErr[3],RMSErr[2],RMSErr[4]);
    else                    printf("Root mean square deviations: P %f, P2 %f, C2 %f \n",RMSErr[0],RMSErr[1],RMSErr[2]);
    
    std::vector<double> Pearson(5,0);   // (qp - q*p) / sqrt(q^2 - (q)^2) sqrt(p^2 - (p)^2)
    for (int i=0;i<Pearson.size();i++) Pearson[i] = (averages[i][2]-(averages[i][0]*averages[i][1])) / (sqrt(averages[i][3]-(averages[i][0]*averages[i][0])) * sqrt(averages[i][4]-(averages[i][1]*averages[i][1])));

    if (computeThreePoints) printf("Pearson correlations: P %f, P2 %f, P3 %f, C2 %f, C3 %f \n",Pearson[0],Pearson[1],Pearson[3],Pearson[2],Pearson[4]);
    else                    printf("Pearson correlations: P %f, P2 %f, C2 %f \n",Pearson[0],Pearson[1],Pearson[2]);
    
    return EXIT_SUCCESS;

}

