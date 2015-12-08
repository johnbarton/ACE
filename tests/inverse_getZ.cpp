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
#include <assert.h>

#include "../src/tools.h"           // Numerical tools
#include "../src/inverse.h"         // Computations of cluster entropies
#include "../src/sparseInverse.h"   // Sparse computations of cluster entropies

int main(int argc, char *argv[]) {

    // Define test spin system (3 sites)
    
    int length  = 3;
    int Jlength = 6;
    Vector J(Jlength), Joff(Jlength), expJ(Jlength), expJoff(Jlength), grad(Jlength), gradoff(Jlength), p(Jlength), poff(Jlength);
    
    int sizes [6] = {1, 2, 3, 2, 3, 6};
    for (int i=0;i<6;i++) {
    
        J[i].resize(sizes[i],0);
        Joff[i].resize(sizes[i],0);
        expJ[i].resize(sizes[i],0);
        expJoff[i].resize(sizes[i],0);
        grad[i].resize(sizes[i],0);
        gradoff[i].resize(sizes[i],0);
        p[i].resize(sizes[i],0);
        poff[i].resize(sizes[i],0);
        
    }

    J[0][0] = 0;
    J[1][0] = 0; J[1][1] = 1;
    J[2][0] = 0; J[2][1] = 1; J[2][2] = 2;
    J[3][0] = 1; J[3][1] = 1;
    J[4][0] =-1; J[4][1] = 0; J[4][2] =-2;
    J[5][0] =-1; J[5][1] = 0; J[5][2] =-2; J[5][3] =-1; J[5][4] = 0; J[5][5] =-2;
    //J[4][0] =-1; J[4][1] =-2; J[4][2] = 0;
    //J[5][0] =-1; J[5][1] =-2; J[5][2] = 0; J[5][3] =-1; J[5][4] =-2; J[5][5] = 0;
    
    Joff[0][0] = 0.5;
    Joff[1][0] = 0;   Joff[1][1] =  0.5;
    Joff[2][0] = 0;   Joff[2][1] =  1;   Joff[2][2] =  0.5;
    Joff[3][0] = 0.5; Joff[3][1] =  0;
    Joff[4][0] = 0.5; Joff[4][1] = -1;   Joff[4][2] = -0.5;
    Joff[5][0] = 0.5; Joff[5][1] = -1;   Joff[5][2] = -2;   Joff[5][3] =  0; Joff[5][4] = 0; Joff[5][5] = 0;
    
    for (int i=0;i<J.size();i++) { for (int a=0;a<J[i].size();a++) {
    
        expJ[i][a]    = exp(J[i][a]);
        expJoff[i][a] = exp(Joff[i][a]);
        
    } }
    
    double Z    = 1;
    double Zoff = 1;
    
    // Note: need as many layers as length when iterating through Z in this way
    
    for (int i=0;i<length;i++) { for (int a=0;a<J[i].size();a++) {
    
        double weight    = expJ[i][a];
        double weightoff = expJoff[i][a];
    
        p[i][a]    += weight;
        Z          += weight;
        
        poff[i][a] += weightoff;
        Zoff       += weightoff;
        
        for (int j=i+1;j<length;j++) { for (int b=0;b<J[j].size();b++) {
        
            int ijx = index(i,j,length);
            int sab = sindex(a,b,J[i].size(),J[j].size());
        
            weight    = expJ[i][a]*expJ[j][b]*expJ[ijx][sab];
            weightoff = expJoff[i][a]*expJoff[j][b]*expJoff[ijx][sab];
            
            p[i][a]        += weight;
            p[j][b]        += weight;
            p[ijx][sab]    += weight;
            Z              += weight;
            
            poff[i][a]     += weightoff;
            poff[j][b]     += weightoff;
            poff[ijx][sab] += weightoff;
            Zoff           += weightoff;
            
            for (int k=j+1;k<length;k++) { for (int c=0;c<J[k].size();c++) {
                
                int ikx = index(i,k,length);
                int jkx = index(j,k,length);
                int sac = sindex(a,c,J[i].size(),J[k].size());
                int sbc = sindex(b,c,J[j].size(),J[k].size());
            
                weight    = expJ[i][a]*expJ[j][b]*expJ[k][c]*expJ[ijx][sab]*expJ[ikx][sac]*expJ[jkx][sbc];
                weightoff = expJoff[i][a]*expJoff[j][b]*expJoff[k][c]*expJoff[ijx][sab]*expJoff[ikx][sac]*expJoff[jkx][sbc];
                
                p[i][a]        += weight;
                p[j][b]        += weight;
                p[k][c]        += weight;
                p[ijx][sab]    += weight;
                p[ikx][sac]    += weight;
                p[jkx][sbc]    += weight;
                Z              += weight;
                
                poff[i][a]     += weightoff;
                poff[j][b]     += weightoff;
                poff[k][c]     += weightoff;
                poff[ijx][sab] += weightoff;
                poff[ikx][sac] += weightoff;
                poff[jkx][sbc] += weightoff;
                Zoff           += weightoff;
                
            } }
            
        } }
    
    } }
    
    for (int i=0;i<p.size();i++) { for (int a=0;a<p[i].size();a++) {
    
        p[i][a]    /= Z;
        poff[i][a] /= Zoff;
        
    } }
    
    for (int i=0;i<p.size();i++) { for (int a=0;a<p[i].size();a++) {
    
        grad[i][a]    = 0;
        gradoff[i][a] = poff[i][a]-p[i][a];
        
    } }

    // Set variables
    
    double tol = 1e-14;
    Vector gradtest(grad);
    std::vector<int> bt;
    std::vector<int> spins;
    
    // Test getZ_gradOnly (not sparse)
    
    double Ztest = getZ_gradOnly(length, expJ, bt, spins, gradtest, length, 1.0);
    for (int i=0;i<J.size();i++) { for (int a=0;a<gradtest[i].size();a++) gradtest[i][a] /= Ztest; }
    bt.clear();
    spins.clear();
    
    if (   fabs(Ztest-Z)>tol) printf("partition function\tgot: %.16e\texpected: %.16e\tdifference: %.16e\n",Ztest,Z,Ztest-Z);
    assert(fabs(Ztest-Z)<tol && "inverse - getZ_gradOnly gave an unexpected result (incorrect partition function)");
    for (int i=0;i<p.size();i++) { for (int a=0;a<p[i].size();a++) {
    
        if (   fabs(p[i][a]-gradtest[i][a])>tol) printf("i_index %d, a_index %d\tgot: %.16e\texpected: %.16e\tdifference: %.16e\n",i,a,gradtest[i][a],p[i][a],gradtest[i][a]-p[i][a]);
        assert(fabs(p[i][a]-gradtest[i][a])<tol && "inverse - getZ_gradOnly gave an unexpected result (incorrect correlation value)");
        
        gradtest[i][a] = 0;
        
    } }
    
    
    Ztest = getZ_gradOnly(length, expJoff, bt, spins, gradtest, length, 1.0);
    for (int i=0;i<J.size();i++) { for (int a=0;a<gradtest[i].size();a++) gradtest[i][a] /= Ztest; }
    bt.clear();
    spins.clear();
    
    if (   fabs(Ztest-Zoff)>tol) printf("partition function\tgot: %.16e\texpected: %.16e\tdifference: %.16e\n",Ztest,Zoff,Ztest-Zoff);
    assert(fabs(Ztest-Zoff)<tol && "inverse - getZ_gradOnly gave an unexpected result (incorrect partition function, off)");
    for (int i=0;i<p.size();i++) { for (int a=0;a<p[i].size();a++) {
        
        if (   fabs(poff[i][a]-gradtest[i][a])>tol) printf("i_index %d, a_index %d\tgot: %.16e\texpected: %.16e\tdifference: %.16e\n",i,a,gradtest[i][a],p[i][a],gradtest[i][a]-p[i][a]);
        assert(fabs(poff[i][a]-gradtest[i][a])<tol && "inverse - getZ_gradOnly gave an unexpected result (incorrect correlation value, off)");
        
        gradtest[i][a]=0;
        
    } }


    // Test getZ (not sparse)
    
    int invLength=0;
    Vector hess;
    std::vector<double> linearHess;
    initialize(J, length, expJ, hess, invLength);
    
    for (int i=0;i<hess.size();i++) { for (int j=0;j<hess[i].size();j++) hess[i][j]=0; }
    Ztest = getZ(length, expJ, bt, spins, gradtest, hess, length, 1.0);
    for (int i=0;i<J.size();i++) { for (int a=0;a<gradtest[i].size();a++) gradtest[i][a] /= Ztest; }
    bt.clear();
    spins.clear();
    
    if (   fabs(Ztest-Z)>tol) printf("partition function\tgot: %.16e\texpected: %.16e\tdifference: %.16e\n",Ztest,Z,Ztest-Z);
    assert(fabs(Ztest-Z)<tol && "inverse - getZ gave an unexpected result (incorrect partition function)");
    for (int i=0;i<p.size();i++) { for (int a=0;a<p[i].size();a++) {
    
        if (   fabs(p[i][a]-gradtest[i][a])>tol) printf("i_index %d, a_index %d\tgot: %.16e\texpected: %.16e\tdifference: %.16e\n",i,a,gradtest[i][a],p[i][a],gradtest[i][a]-p[i][a]);
        assert(fabs(p[i][a]-gradtest[i][a])<tol && "inverse - getZ gave an unexpected result (incorrect correlation value)");
        
        gradtest[i][a] = 0;
        
    } }
    
    
    for (int i=0;i<hess.size();i++) { for (int j=0;j<hess[i].size();j++) hess[i][j]=0; }
    Ztest = getZ(length, expJoff, bt, spins, gradtest, hess, length, 1.0);
    for (int i=0;i<J.size();i++) { for (int a=0;a<gradtest[i].size();a++) gradtest[i][a] /= Ztest; }
    bt.clear();
    spins.clear();
    
    if (   fabs(Ztest-Zoff)>tol) printf("partition function\tgot: %.16e\texpected: %.16e\tdifference: %.16e\n",Ztest,Zoff,Ztest-Zoff);
    assert(fabs(Ztest-Zoff)<tol && "inverse - getZ gave an unexpected result (incorrect partition function, off)");
    for (int i=0;i<p.size();i++) { for (int a=0;a<p[i].size();a++) {
    
        if (   fabs(poff[i][a]-gradtest[i][a])>tol) printf("i_index %d, a_index %d\tgot: %.16e\texpected: %.16e\tdifference: %.16e\n",i,a,gradtest[i][a],p[i][a],gradtest[i][a]-p[i][a]);
        assert(fabs(poff[i][a]-gradtest[i][a])<tol && "inverse - getZ gave an unexpected result (incorrect correlation value, off)");
        
        gradtest[i][a] = 0;
        
    } }

    
    
    
    // Set variables for sparse routines

    IntVector interactingSpins(length,std::vector<int>());
    IntVector independentSpins(length,std::vector<int>());
    IntVector interactingSpinsoff(length,std::vector<int>());
    IntVector independentSpinsoff(length,std::vector<int>());
    
    interactingSpins[0].push_back(0);
    interactingSpins[1].push_back(0); interactingSpins[1].push_back(1);
    interactingSpins[2].push_back(0); independentSpins[2].push_back(1); interactingSpins[2].push_back(2);
    //interactingSpins[2].push_back(0); interactingSpins[2].push_back(1); independentSpins[2].push_back(2);
    
    interactingSpinsoff[0].push_back(0);
    interactingSpinsoff[1].push_back(0); interactingSpinsoff[1].push_back(1);
    interactingSpinsoff[2].push_back(0); interactingSpinsoff[2].push_back(1); interactingSpinsoff[2].push_back(2);
    
    
    // Test getZ_gradOnly (sparse)
    
    int btarr[length];
    int btC[length];
    int spinsarr[length];
    for (int i=0;i<length;i++) { btarr[i]=length; btC[i]=i; spinsarr[i]=0; }
    
    double stem=1;
    double zero[length];
    for (int i=0;i<length;i++) {
    
        zero[i] = 1;
        for (int j=0;j<independentSpins[i].size();j++) zero[i] += expJ[i][independentSpins[i][j]];
        stem *= zero[i];
        
    }
    
    /* optimizeS code */
    
            for (int i=0;i<length;++i) { for (int j=0;j<independentSpins[btC[i]].size();j++) {
            
                double ijStem = stem * expJ[btC[i]][independentSpins[btC[i]][j]] / zero[btC[i]];
                
                gradtest[btC[i]][independentSpins[btC[i]][j]] += ijStem;
                
                int off = offset(btC[i],length);
                
                for (int k=i+1;k<length;++k) { for (int l=0;l<independentSpins[btC[k]].size();l++) {
                
                    int sik = sindex(independentSpins[btC[i]][j],independentSpins[btC[k]][l],expJ[btC[i]].size(),expJ[btC[k]].size());
                    
                    gradtest[off + btC[k]][sik] += ijStem * expJ[btC[k]][independentSpins[btC[k]][l]] / zero[btC[k]];
                    
                } }
                
            } }
    
            Ztest = stem;
    
            for (int i=0;i<length;i++) {
            
                eraseInPlace(btC,length,i);
                
                for (int p=0;p<interactingSpins[i].size();p++) {
                    
                    double prodJ = expJ[i][interactingSpins[i][p]];
                    
                    btarr[0]=i;
                    spinsarr[i]=interactingSpins[i][p];
                    
                    Ztest += getZ_gradOnly(1, expJ, btarr, btC, spinsarr, gradtest, length, prodJ * stem / zero[i], zero, interactingSpins, independentSpins);
                
                }
                
                insertInPlace(btC,length,i);
                
            }
            
            for (int i=0;i<J.size();i++) { for (int j=0;j<gradtest[i].size();j++) gradtest[i][j] /= Ztest; }
    
    /* end optimizeS code */
    
    if (   fabs(Ztest-Z)>tol) printf("partition function\tgot: %.16e\texpected: %.16e\tdifference: %.16e\n",Ztest,Z,Ztest-Z);
    assert(fabs(Ztest-Z)<tol && "inverse - getZ_gradOnly, sparse gave an unexpected result (incorrect partition function)");
    for (int i=0;i<p.size();i++) { for (int a=0;a<p[i].size();a++) {
    
        if (   fabs(p[i][a]-gradtest[i][a])>tol) printf("i_index %d, a_index %d\tgot: %.16e\texpected: %.16e\tdifference: %.16e\n",i,a,gradtest[i][a],p[i][a],gradtest[i][a]-p[i][a]);
        assert(fabs(p[i][a]-gradtest[i][a])<tol && "inverse - getZ_gradOnly, sparse gave an unexpected result (incorrect correlation value)");
        
        gradtest[i][a] = 0;
        
    } }

    
    // Test getZ (sparse) NOTE: HESSIAN IS IGNORED
    
    for (int i=0;i<length;i++) { btarr[i]=length; btC[i]=i; spinsarr[i]=0; }
    stem=1;
    for (int i=0;i<length;i++) {
    
        zero[i] = 1;
        for (int j=0;j<independentSpins[i].size();j++) zero[i] += expJ[i][independentSpins[i][j]];
        stem *= zero[i];
        
    }
    
    /* optimizeS code */
    
            for (int i=0;i<length;++i) { for (int j=0;j<independentSpins[btC[i]].size();j++) {
            
                double ijStem = stem * expJ[btC[i]][independentSpins[btC[i]][j]] / zero[btC[i]];
                
                gradtest[btC[i]][independentSpins[btC[i]][j]] += ijStem;
                
                int off = offset(btC[i],length);
                
                for (int k=i+1;k<length;++k) { for (int l=0;l<independentSpins[btC[k]].size();l++) {
                
                    int sik = sindex(independentSpins[btC[i]][j],independentSpins[btC[k]][l],expJ[btC[i]].size(),expJ[btC[k]].size());
                    
                    gradtest[off + btC[k]][sik] += ijStem * expJ[btC[k]][independentSpins[btC[k]][l]] / zero[btC[k]];
                    
                } }
                
            } }
    
            Ztest = stem;
    
            for (int i=0;i<length;i++) {
    
                eraseInPlace(btC,length,i);
                
                for (int p=0;p<interactingSpins[i].size();p++) {
                    
                    double prodJ = expJ[i][interactingSpins[i][p]];
                    
                    btarr[0]=i;
                    spinsarr[i]=interactingSpins[i][p];
                    
                    Ztest += getZ(1, expJ, btarr, btC, spinsarr, gradtest, hess, length, prodJ * stem / zero[i], zero, interactingSpins, independentSpins);
                
                }
                
                insertInPlace(btC,length,i);
                
            }
            
            for (int i=0;i<J.size();i++) { for (int j=0;j<gradtest[i].size();j++) gradtest[i][j] /= Ztest; }
    
    /* end optimizeS code */
    
    if (   fabs(Ztest-Z)>tol) printf("partition function\tgot: %.16e\texpected: %.16e\tdifference: %.16e\n",Ztest,Z,Ztest-Z);
    assert(fabs(Ztest-Z)<tol && "inverse - getZ, sparse gave an unexpected result (incorrect partition function)");
    for (int i=0;i<p.size();i++) { for (int a=0;a<p[i].size();a++) {
    
        if (   fabs(p[i][a]-gradtest[i][a])>tol) printf("i_index %d, a_index %d\tgot: %.16e\texpected: %.16e\tdifference: %.16e\n",i,a,gradtest[i][a],p[i][a],gradtest[i][a]-p[i][a]);
        assert(fabs(p[i][a]-gradtest[i][a])<tol && "inverse - getZ, sparse gave an unexpected result (incorrect correlation value)");
        
        gradtest[i][a] = 0;
        
    } }
    
    
    
    
    // Test sparse routines with dense couplings
    
    // Test getZ_gradOnly (sparse)
    
    for (int i=0;i<length;i++) { btarr[i]=length; btC[i]=i; spinsarr[i]=0; }
    
    stem=1;
    for (int i=0;i<length;i++) {
    
        zero[i] = 1;
        for (int j=0;j<independentSpinsoff[i].size();j++) zero[i] += expJ[i][independentSpinsoff[i][j]];
        stem *= zero[i];
        
    }
    
    /* optimizeS code */
    
            for (int i=0;i<length;++i) { for (int j=0;j<independentSpinsoff[btC[i]].size();j++) {
            
                double ijStem = stem * expJ[btC[i]][independentSpinsoff[btC[i]][j]] / zero[btC[i]];
                
                gradtest[btC[i]][independentSpinsoff[btC[i]][j]] += ijStem;
                
                int off = offset(btC[i],length);
                
                for (int k=i+1;k<length;++k) { for (int l=0;l<independentSpinsoff[btC[k]].size();l++) {
                
                    int sik = sindex(independentSpinsoff[btC[i]][j],independentSpinsoff[btC[k]][l],expJ[btC[i]].size(),expJ[btC[k]].size());
                    
                    gradtest[off + btC[k]][sik] += ijStem * expJ[btC[k]][independentSpinsoff[btC[k]][l]] / zero[btC[k]];
                    
                } }
                
            } }
    
            Ztest = stem;
    
            for (int i=0;i<length;i++) {
            
                eraseInPlace(btC,length,i);
                
                for (int p=0;p<interactingSpinsoff[i].size();p++) {
                    
                    double prodJ = expJ[i][interactingSpinsoff[i][p]];
                    
                    btarr[0]=i;
                    spinsarr[i]=interactingSpinsoff[i][p];
                    
                    Ztest += getZ_gradOnly(1, expJ, btarr, btC, spinsarr, gradtest, length, prodJ * stem / zero[i], zero, interactingSpinsoff, independentSpinsoff);
                
                }
                
                insertInPlace(btC,length,i);
                
            }
            
            for (int i=0;i<J.size();i++) { for (int j=0;j<gradtest[i].size();j++) gradtest[i][j] /= Ztest; }
    
    /* end optimizeS code */
    
    if (   fabs(Ztest-Z)>tol) printf("partition function\tgot: %.16e\texpected: %.16e\tdifference: %.16e\n",Ztest,Z,Ztest-Z);
    assert(fabs(Ztest-Z)<tol && "inverse - getZ_gradOnly, sparse gave an unexpected result (incorrect partition function)");
    for (int i=0;i<p.size();i++) { for (int a=0;a<p[i].size();a++) {
    
        if (   fabs(p[i][a]-gradtest[i][a])>tol) printf("i_index %d, a_index %d\tgot: %.16e\texpected: %.16e\tdifference: %.16e\n",i,a,gradtest[i][a],p[i][a],gradtest[i][a]-p[i][a]);
        assert(fabs(p[i][a]-gradtest[i][a])<tol && "inverse - getZ_gradOnly, sparse gave an unexpected result (incorrect correlation value)");
        
        gradtest[i][a] = 0;
        
    } }
    
    
    // Test getZ (sparse) NOTE: HESSIAN IS IGNORED
    
    for (int i=0;i<length;i++) { btarr[i]=length; btC[i]=i; spinsarr[i]=0; }
    stem=1;
    for (int i=0;i<length;i++) {
    
        zero[i] = 1;
        for (int j=0;j<independentSpinsoff[i].size();j++) zero[i] += expJ[i][independentSpinsoff[i][j]];
        stem *= zero[i];
        
    }
    
    /* optimizeS code */
    
            for (int i=0;i<length;++i) { for (int j=0;j<independentSpinsoff[btC[i]].size();j++) {
            
                double ijStem = stem * expJ[btC[i]][independentSpinsoff[btC[i]][j]] / zero[btC[i]];
                
                gradtest[btC[i]][independentSpinsoff[btC[i]][j]] += ijStem;
                
                int off = offset(btC[i],length);
                
                for (int k=i+1;k<length;++k) { for (int l=0;l<independentSpinsoff[btC[k]].size();l++) {
                
                    int sik = sindex(independentSpinsoff[btC[i]][j],independentSpinsoff[btC[k]][l],expJ[btC[i]].size(),expJ[btC[k]].size());
                    
                    gradtest[off + btC[k]][sik] += ijStem * expJ[btC[k]][independentSpinsoff[btC[k]][l]] / zero[btC[k]];
                    
                } }
                
            } }
    
            Ztest = stem;
    
            for (int i=0;i<length;i++) {
    
                eraseInPlace(btC,length,i);
                
                for (int p=0;p<interactingSpinsoff[i].size();p++) {
                    
                    double prodJ = expJ[i][interactingSpinsoff[i][p]];
                    
                    btarr[0]=i;
                    spinsarr[i]=interactingSpinsoff[i][p];
                    
                    Ztest += getZ(1, expJ, btarr, btC, spinsarr, gradtest, hess, length, prodJ * stem / zero[i], zero, interactingSpinsoff, independentSpinsoff);
                
                }
                
                insertInPlace(btC,length,i);
                
            }
            
            for (int i=0;i<J.size();i++) { for (int j=0;j<gradtest[i].size();j++) gradtest[i][j] /= Ztest; }
    
    /* end optimizeS code */
    
    if (   fabs(Ztest-Z)>tol) printf("partition function\tgot: %.16e\texpected: %.16e\tdifference: %.16e\n",Ztest,Z,Ztest-Z);
    assert(fabs(Ztest-Z)<tol && "inverse - getZ, sparse gave an unexpected result (incorrect partition function)");
    for (int i=0;i<p.size();i++) { for (int a=0;a<p[i].size();a++) {
    
        if (   fabs(p[i][a]-gradtest[i][a])>tol) printf("i_index %d, a_index %d\tgot: %.16e\texpected: %.16e\tdifference: %.16e\n",i,a,gradtest[i][a],p[i][a],gradtest[i][a]-p[i][a]);
        assert(fabs(p[i][a]-gradtest[i][a])<tol && "inverse - getZ, sparse gave an unexpected result (incorrect correlation value)");
        
        gradtest[i][a] = 0;
        
    } }
    
    
	return 0;

}

