#include <vector>
#include <set>

#include "tools.h"          // Numerical tools
#include "sparseInverse.h"  // Computations of cluster entropies


//#define EPSG 1.0e-10
//#define MAXANGLE 0.998        
//#define MAXSTEPSIZE 5.0
//#define MINSTEPSIZE 1.0e-10
//#define NEWTONGRADSIZE  3.0e-3
//#define NEWTONGRADRATIO 1.0e-7
//#define ARMIJOC 1.0e-4
//#define ARMIJOP 0.5

#define EPSG 1.0e-9
#define MAXANGLE 0.998        
#define MAXSTEPSIZE 5.0
#define MINSTEPSIZE 1.0e-11
#define NEWTONGRADSIZE  3.0e-3
#define NEWTONGRADRATIO 1.0e-8
#define ARMIJOC 1.0e-5
#define ARMIJOP 0.5
#define STEMCUT 1.0e-14





/*****************************************
       
   P A R T I T I O N    F U N C T I O N
 
*****************************************/





// Recursive computation of the partition function

double getZ(int L, const Vector &expJ, int bt[], int btC[], int spins[], Vector &grad, Vector &hess, int length, double stem, double zero[], const IntVector &interactingSpins, const IntVector &independentSpins) {
        
    // Compute correlations for independent spins
    
    for (int i=0;i<length-L;++i) { for (int j=0;j<independentSpins[btC[i]].size();j++) {
    
        double ijStem = stem * expJ[btC[i]][independentSpins[btC[i]][j]] / zero[btC[i]];
        
        grad[btC[i]][independentSpins[btC[i]][j]] += ijStem;
        hess[hindex(btC[i],btC[i],expJ.size())][hindex(independentSpins[btC[i]][j],independentSpins[btC[i]][j],expJ[btC[i]].size())] += ijStem;
            
        int off = offset(btC[i],length);
        
        for (int k=i+1;k<length-L;++k) { for (int l=0;l<independentSpins[btC[k]].size();l++) {
        
            int sik = sindex(independentSpins[btC[i]][j],independentSpins[btC[k]][l],expJ[btC[i]].size(),expJ[btC[k]].size());
        
            grad[off + btC[k]][sik]                      += ijStem * expJ[btC[k]][independentSpins[btC[k]][l]] / zero[btC[k]];
            hess[hindex(btC[i],btC[k],expJ.size())][sik] += ijStem * expJ[btC[k]][independentSpins[btC[k]][l]] / zero[btC[k]];
            
        } }
        
    } }
    
    // Compute correlations for interacting spins
    
    for (int i=0;i<L;++i) {
            
        // One-body
        
        grad[bt[i]][spins[bt[i]]] += stem;
        hess[hindex(bt[i],bt[i],expJ.size())][hindex(spins[bt[i]],spins[bt[i]],expJ[bt[i]].size())] += stem;
        
        // Two-Body (interacting)
        
        int offi = offset(bt[i],length);
        
        for (int j=i+1;j<L;++j) {
        
            int sij = sindex(spins[bt[i]],spins[bt[j]],expJ[bt[i]].size(),expJ[bt[j]].size());
        
            grad[offi + bt[j]][sij] += stem;
            
            hess[hindex(bt[i],bt[j],expJ.size())][sij] += stem;
            hess[hindex(bt[i],offi+bt[j],expJ.size())][sindex(spins[bt[i]],sij,expJ[bt[i]].size(),expJ[offi+bt[j]].size())] += stem;
            hess[hindex(bt[j],offi+bt[j],expJ.size())][sindex(spins[bt[j]],sij,expJ[bt[j]].size(),expJ[offi+bt[j]].size())] += stem;
            hess[hindex(offi+bt[j],offi+bt[j],expJ.size())][hindex(sij,sij,expJ[offi+bt[j]].size())] += stem;
            
            // Three-body (interacting)
            
            int offj = offset(bt[j],length);
            
            for (int k=j+1;k<L;++k) {
            
                int sik = sindex(spins[bt[i]],spins[bt[k]],expJ[bt[i]].size(),expJ[bt[k]].size());
                int sjk = sindex(spins[bt[j]],spins[bt[k]],expJ[bt[j]].size(),expJ[bt[k]].size());
                
                hess[hindex(bt[i],offj+bt[k],expJ.size())][sindex(spins[bt[i]],sjk,expJ[bt[i]].size(),expJ[offj+bt[k]].size())] += stem;
                hess[hindex(bt[j],offi+bt[k],expJ.size())][sindex(spins[bt[j]],sik,expJ[bt[j]].size(),expJ[offi+bt[k]].size())] += stem;
                hess[hindex(bt[k],offi+bt[j],expJ.size())][sindex(spins[bt[k]],sij,expJ[bt[k]].size(),expJ[offi+bt[j]].size())] += stem;
                    
                hess[hindex(offi+bt[j],offi+bt[k],expJ.size())][sindex(sij,sik,expJ[offi+bt[j]].size(),expJ[offi+bt[k]].size())] += stem;
                hess[hindex(offi+bt[j],offj+bt[k],expJ.size())][sindex(sij,sjk,expJ[offi+bt[j]].size(),expJ[offj+bt[k]].size())] += stem;
                hess[hindex(offi+bt[k],offj+bt[k],expJ.size())][sindex(sik,sjk,expJ[offi+bt[k]].size(),expJ[offj+bt[k]].size())] += stem;
                
                // Four-body (interacting)
                
                int offk=offset(bt[k],length);
                    
                for (int l=k+1;l<L;++l) {
                    
                    int sil = sindex(spins[bt[i]],spins[bt[l]],expJ[bt[i]].size(),expJ[bt[l]].size());
                    int sjl = sindex(spins[bt[j]],spins[bt[l]],expJ[bt[j]].size(),expJ[bt[l]].size());
                    int skl = sindex(spins[bt[k]],spins[bt[l]],expJ[bt[k]].size(),expJ[bt[l]].size());
                        
                    hess[hindex(offi+bt[j],offk+bt[l],expJ.size())][sindex(sij,skl,expJ[offi+bt[j]].size(),expJ[offk+bt[l]].size())] += stem;
                    hess[hindex(offi+bt[k],offj+bt[l],expJ.size())][sindex(sik,sjl,expJ[offi+bt[k]].size(),expJ[offj+bt[l]].size())] += stem;
                    hess[hindex(offi+bt[l],offj+bt[k],expJ.size())][sindex(sil,sjk,expJ[offi+bt[l]].size(),expJ[offj+bt[k]].size())] += stem;
                    
                }
            
            }
            
            // Three-body (independent)
            
            for (int k=0;k<length-L;++k) { for (int l=0;l<independentSpins[btC[k]].size();l++) {
                    
                hess[hindex(btC[k],offi+bt[j],expJ.size())][sindex(independentSpins[btC[k]][l],sij,expJ[btC[k]].size(),expJ[offi+bt[j]].size())] += stem * expJ[btC[k]][independentSpins[btC[k]][l]] / zero[btC[k]];
                    
            } }
            
        }
        
        // Two-body (independent)
        
        int j;
        for (j=0;j<length-L && btC[j]<bt[i];++j) { for (int k=0;k<independentSpins[btC[j]].size();k++) {
        
            int sji = sindex(independentSpins[btC[j]][k],spins[bt[i]],expJ[btC[j]].size(),expJ[bt[i]].size());
            
            double jkStem = stem * expJ[btC[j]][independentSpins[btC[j]][k]] / zero[btC[j]];
                
            grad[index(btC[j],bt[i],length)][sji]       += jkStem;
            hess[hindex(btC[j],bt[i],expJ.size())][sji] += jkStem;
                    
        } }
        
        for (;j<length-L;++j) { for (int k=0;k<independentSpins[btC[j]].size();k++) {
        
            int sij = sindex(spins[bt[i]],independentSpins[btC[j]][k],expJ[bt[i]].size(),expJ[btC[j]].size());
            
            double jkStem = stem * expJ[btC[j]][independentSpins[btC[j]][k]] / zero[btC[j]];
                
            grad[index(bt[i],btC[j],length)][sij]       += jkStem;
            hess[hindex(bt[i],btC[j],expJ.size())][sij] += jkStem;
                    
        } }
    
    }
    
    // Truncation condition
    
    //if (stem<STEMCUT) return stem;
    
    // Branch further down the tree
    
    double temp = stem;
    
    for (int i=bt[L-1]+1;i<length;++i) {
        
        eraseInPlace(btC,length-L,i);
        
        for (int p=0;p<interactingSpins[i].size();p++) {
            
            double prodJ = expJ[i][interactingSpins[i][p]];
            
            for (int j=0;j<L;++j) prodJ *= expJ[index(bt[j],i,length)][sindex(spins[bt[j]],interactingSpins[i][p],expJ[bt[j]].size(),expJ[i].size())];
            
            bt[L]=i;
            spins[i]=interactingSpins[i][p];
            
            temp += getZ(L+1, expJ, bt, btC, spins, grad, hess, length, prodJ * stem / zero[i], zero, interactingSpins, independentSpins);
        
        }
        
        insertInPlace(btC,length-L,i);
    
    }
    
    // Pass up the tree
    
    return temp;

}


// Recursive computation of the cross entropy

void optimizeS(const Vector &J, const Vector &expJ, double &S, Vector &grad, Vector &hess, const Vector &p, const IntVector &interactingSpins, const IntVector &independentSpins) {

    // Set variables
    
    int length=(int) sizetolength(J.size());

    for (int i=0;i<grad.size();i++) { for (int j=0;j<grad[i].size();j++) grad[i][j]=0; }
    for (int i=0;i<hess.size();i++) { for (int j=0;j<hess[i].size();j++) hess[i][j]=0; }
    
    int bt[length];
    int btC[length];
    int spins[length];
    for (int i=0;i<length;i++) { bt[i]=length; btC[i]=i; spins[i]=0; }
    
    double stem=1;
    double zero[length];
    for (int i=0;i<length;i++) {
    
        zero[i] = 1;
        for (int j=0;j<independentSpins[i].size();j++) zero[i] += expJ[i][independentSpins[i][j]];
        stem *= zero[i];
        
    }
    
    // Run first iteration over all zero states
    
    for (int i=0;i<length;++i) { for (int j=0;j<independentSpins[btC[i]].size();j++) {
    
        double ijStem = stem * expJ[btC[i]][independentSpins[btC[i]][j]] / zero[btC[i]];
        
        grad[btC[i]][independentSpins[btC[i]][j]] += ijStem;
        hess[hindex(btC[i],btC[i],expJ.size())][hindex(independentSpins[btC[i]][j],independentSpins[btC[i]][j],expJ[btC[i]].size())] += ijStem;
            
        int off = offset(btC[i],length);
        
        for (int k=i+1;k<length;++k) { for (int l=0;l<independentSpins[btC[k]].size();l++) {
        
            int sik = sindex(independentSpins[btC[i]][j],independentSpins[btC[k]][l],expJ[btC[i]].size(),expJ[btC[k]].size());
            
            grad[off + btC[k]][sik]                      += ijStem * expJ[btC[k]][independentSpins[btC[k]][l]] / zero[btC[k]];
            hess[hindex(btC[i],btC[k],expJ.size())][sik] += ijStem * expJ[btC[k]][independentSpins[btC[k]][l]] / zero[btC[k]];
            
        } }
        
    } }
    
    // Get contributions to partition function

    double Z = stem;
    
    for (int i=0;i<length;i++) {
    
        eraseInPlace(btC,length,i);
        
        for (int p=0;p<interactingSpins[i].size();p++) {
            
            double prodJ = expJ[i][interactingSpins[i][p]];
            
            bt[0]=i;
            spins[i]=interactingSpins[i][p];
            
            Z += getZ(1, expJ, bt, btC, spins, grad, hess, length, prodJ * stem / zero[i], zero, interactingSpins, independentSpins);
        
        }
        
        insertInPlace(btC,length,i);
        
    }
    
    // Get correlations
    
    for (int i=0;i<J.size();i++) { for (int j=0;j<grad[i].size();j++) grad[i][j] /= Z; }
    for (int i=0;i<J.size();i++) { for (int j=0;j<grad[i].size();j++) {
        
        for (int l=j;l<grad[i].size();l++) {
        
            hess[hindex(i,i,J.size())][hindex(j,l,grad[i].size())] /= Z;
            hess[hindex(i,i,J.size())][hindex(j,l,grad[i].size())] -= grad[i][j] * grad[i][l];
            
        }
    
        for (int k=i+1;k<J.size();k++) { for (int l=0;l<grad[k].size();l++) {
        
            hess[hindex(i,k,J.size())][sindex(j,l,grad[i].size(),grad[k].size())] /= Z;
            hess[hindex(i,k,J.size())][sindex(j,l,grad[i].size(),grad[k].size())] -= grad[i][j] * grad[k][l];
        
        } }
        
    } }
    
    // Calculate p * J
	
	double p_dot_J=0;
    
	for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) {
        
        p_dot_J    += J[i][j] * p[i][j];
        grad[i][j] -= p[i][j];
        
    } }
	
	// Return the result
	
	S = (log(Z) - p_dot_J);

}


// Recursive computation of the partition function (gradient only)

double getZ_gradOnly(int L, const Vector &expJ, int bt[], int btC[], int spins[], Vector &grad, int length, double stem, double zero[], const IntVector &interactingSpins, const IntVector &independentSpins) {
        
    // Compute correlations for independent spins
    
    for (int i=0;i<length-L;++i) { for (int j=0;j<independentSpins[btC[i]].size();j++) {
    
        double ijStem = stem * expJ[btC[i]][independentSpins[btC[i]][j]] / zero[btC[i]];
        
        grad[btC[i]][independentSpins[btC[i]][j]] += ijStem;
        
        int off = offset(btC[i],length);
        
        for (int k=i+1;k<length-L;++k) { for (int l=0;l<independentSpins[btC[k]].size();l++) {
        
            int sik = sindex(independentSpins[btC[i]][j],independentSpins[btC[k]][l],expJ[btC[i]].size(),expJ[btC[k]].size());
        
            grad[off + btC[k]][sik] += ijStem * expJ[btC[k]][independentSpins[btC[k]][l]] / zero[btC[k]];
            
        } }
        
    } }
    
    // Compute correlations for interacting spins
    
    for (int i=0;i<L;++i) {
            
        // One-body
        
        grad[bt[i]][spins[bt[i]]] += stem;
        
        // Two-Body (interacting)
        
        int offi = offset(bt[i],length);
        
        for (int j=i+1;j<L;++j) {
        
            int sij = sindex(spins[bt[i]],spins[bt[j]],expJ[bt[i]].size(),expJ[bt[j]].size());
        
            grad[offi + bt[j]][sij] += stem;
        
        }
        
        // Two-body (independent)
        
        int j;
        for (j=0;j<length-L && btC[j]<bt[i];++j) { for (int k=0;k<independentSpins[btC[j]].size();k++) {
        
            int sji = sindex(independentSpins[btC[j]][k],spins[bt[i]],expJ[btC[j]].size(),expJ[bt[i]].size());
        
            double jkStem = stem * expJ[btC[j]][independentSpins[btC[j]][k]] / zero[btC[j]];
                
            grad[index(btC[j],bt[i],length)][sji] += jkStem;
            
        } }
        
        for (;j<length-L;++j) { for (int k=0;k<independentSpins[btC[j]].size();k++) {
        
            int sij = sindex(spins[bt[i]],independentSpins[btC[j]][k],expJ[bt[i]].size(),expJ[btC[j]].size());
            
            double jkStem = stem * expJ[btC[j]][independentSpins[btC[j]][k]] / zero[btC[j]];
                
            grad[index(bt[i],btC[j],length)][sij] += jkStem;
            
        } }
    
    }
    
    // Truncation condition
    
    //if (stem<STEMCUT) return stem;
    
    // Branch further down the tree
    
    double temp = stem;
    
    for (int i=bt[L-1]+1; i<length; ++i) {
        
        eraseInPlace(btC,length-L,i);
        
        for (int p=0;p<interactingSpins[i].size();p++) {
            
            double prodJ = expJ[i][interactingSpins[i][p]];
            
            for (int j=0;j<L;++j) prodJ *= expJ[index(bt[j],i,length)][sindex(spins[bt[j]],interactingSpins[i][p],expJ[bt[j]].size(),expJ[i].size())];
            
            bt[L]=i;
            spins[i]=interactingSpins[i][p];
            
            temp += getZ_gradOnly(L+1, expJ, bt, btC, spins, grad, length, prodJ * stem / zero[i], zero, interactingSpins, independentSpins);
            
            bt[L]=0;
            spins[i]=0;
        
        }
        
        insertInPlace(btC,length-L,i);
    
    }
    
    // Pass up the tree
    
    return temp;

}


// Recursive computation of the cross entropy (gradient only)

void optimizeS_gradOnly(const Vector &J, const Vector &expJ, double &S, Vector &grad, const Vector &p, const IntVector &interactingSpins, const IntVector &independentSpins) {

    // Set variables
    
    int length=(int) sizetolength(J.size());

    for (int i=0;i<grad.size();i++) { for (int j=0;j<grad[i].size();j++) grad[i][j]=0; }
    
    // Get partition function and correlations
    
    int bt[length];
    int btC[length];
    int spins[length];
    for (int i=0;i<length;i++) { bt[i]=length; btC[i]=i; spins[i]=0; }
    
    double stem=1;
    double zero[length];
    for (int i=0;i<length;i++) {
    
        zero[i] = 1;
        for (int j=0;j<independentSpins[i].size();j++) zero[i] += expJ[i][independentSpins[i][j]];
        stem *= zero[i];
        
    }
    
    // Run first iteration over all zero states manually
    
    for (int i=0;i<length;++i) { for (int j=0;j<independentSpins[btC[i]].size();j++) {
    
        double ijStem = stem * expJ[btC[i]][independentSpins[btC[i]][j]] / zero[btC[i]];
        
        grad[btC[i]][independentSpins[btC[i]][j]] += ijStem;
        
        int off = offset(btC[i],length);
        
        for (int k=i+1;k<length;++k) { for (int l=0;l<independentSpins[btC[k]].size();l++) {
        
            int sik = sindex(independentSpins[btC[i]][j],independentSpins[btC[k]][l],expJ[btC[i]].size(),expJ[btC[k]].size());
            
            grad[off + btC[k]][sik] += ijStem * expJ[btC[k]][independentSpins[btC[k]][l]] / zero[btC[k]];
            
        } }
        
    } }
    
    // Get contributions to partition function

    double Z = stem;
    
    for (int i=0;i<length;i++) {
    
        eraseInPlace(btC,length,i);
        
        for (int p=0;p<interactingSpins[i].size();p++) {
            
            double prodJ = expJ[i][interactingSpins[i][p]];
            
            bt[0]=i;
            spins[i]=interactingSpins[i][p];
            
            Z += getZ_gradOnly(1, expJ, bt, btC, spins, grad, length, prodJ * stem / zero[i], zero, interactingSpins, independentSpins);
        
        }
        
        insertInPlace(btC,length,i);
        
    }
    
    for (int i=0;i<J.size();i++) { for (int j=0;j<grad[i].size();j++) grad[i][j] /= Z; }
    
    // Calculate p * J
	
	double p_dot_J=0;
    
	for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) {
        
        p_dot_J    += J[i][j] * p[i][j];
        grad[i][j] -= p[i][j];
        
    } }
	
	// Return the result
	
	S = (log(Z) - p_dot_J);

}





/*****************************************
       
  A U X I L I A R Y    F U N C T I O N S
 
*****************************************/





// Initialize variables for the optimization routine

void initialize(const Vector &J, int length, Vector &expJ, Vector &hess, int &invLength, const IntVector &workingSet, IntVector &interactingSpins, IntVector &independentSpins) {

    invLength=0;

    for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) expJ[i][j] = exp(J[i][j]); }
    
    hess.resize( (J.size() * (J.size() + 1)) / 2 );
    
    for (int i=0;i<length;i++) {
    
        hess[hindex(i,i,J.size())].resize( (J[i].size() * (J[i].size() + 1)) / 2, 0);                      // Resize field-field diagonal
        
        for (int j=i+1;j<length;j++)      hess[hindex(i,j,J.size())].resize(J[i].size() * J[j].size(), 0); // Resize field-field off-diagonal
        for (int j=length;j<J.size();j++) hess[hindex(i,j,J.size())].resize(J[i].size() * J[j].size(), 0); // Resize field-coupling
        
        invLength += (int) workingSet[i].size();
        
    }
    for (int i=length;i<J.size();i++) {
    
        hess[hindex(i,i,J.size())].resize( (J[i].size() * (J[i].size() + 1)) / 2, 0);                      // Resize coupling-coupling diagonal
        
        for (int j=i+1;j<J.size();j++)    hess[hindex(i,j,J.size())].resize(J[i].size() * J[j].size(), 0); // Resize coupling-coupling off-diagonal
        
        invLength += (int) workingSet[i].size();
        
    }
    
    interactingSpins.resize(length,std::vector<int>());
    independentSpins.resize(length,std::vector<int>());
    
    for (int site1=0;site1<length;site1++) { for (int site2=site1+1;site2<length;site2++) {
    
        for (int x=0;x<workingSet[index(site1,site2,length)].size();x++) {
        
            int state2=workingSet[index(site1,site2,length)][x]%J[site2].size();
            int state1=(workingSet[index(site1,site2,length)][x]-state2)/J[site2].size();
            
            insertInPlace(interactingSpins[site1],state1);
            insertInPlace(interactingSpins[site2],state2);
        
        }
        
    } }
    
    for (int i=0;i<length;i++) {
    
        int count=0;
        
        for (int j=0;j<J[i].size();j++) {
    
            if (count==interactingSpins[i].size())  independentSpins[i].push_back(j);
            else if (interactingSpins[i][count]==j) count++;
            else                                    independentSpins[i].push_back(j);
            
        }
        
    }

}


// Choose to take a step using Newton method or gradient descent (sparse version)

bool useDescent(const Vector &grad, const IntVector &workingSet) {
    
    double gradMax=0;
    double gradMin=1;
        
    for (int i=0;i<workingSet.size();i++) { for (int j=0;j<workingSet[i].size();j++) {
            
        if (fabs(grad[i][workingSet[i][j]])>gradMax) gradMax=fabs(grad[i][workingSet[i][j]]);
        if (fabs(grad[i][workingSet[i][j]])<gradMin) gradMin=fabs(grad[i][workingSet[i][j]]);
            
    } }
        
    if (gradMin<EPSG) gradMin=EPSG;
    
    return ( gradMax>NEWTONGRADSIZE || (gradMin/gradMax)<NEWTONGRADRATIO );

}


// Compute step size using gradient descent (sparse version)

void computeDescentStep(const Vector &grad, int length, Vector &step, const IntVector &workingSet) {

    double totalGrad=0;
    if (length==1) totalGrad=1.0;
        
    for (int i=0;i<workingSet.size();i++) { for (int j=0;j<workingSet[i].size();j++) { totalGrad                += fabs(grad[i][workingSet[i][j]]);      } }
    for (int i=0;i<workingSet.size();i++) { for (int j=0;j<workingSet[i].size();j++) { step[i][workingSet[i][j]] = -grad[i][workingSet[i][j]]/totalGrad; } }
    
}


// Compute step size using Newton's method (sparse version)

void computeNewtonStep(const Vector &grad, Vector &hess, std::vector<double> &linearHess, std::vector<double> &inverseC, std::vector<double> &inverseN, std::vector<double> &inverseD, Vector &step, const IntVector &workingSet) {

    modifiedCholeskyInverse_sparse(hess,linearHess,inverseC,inverseN,inverseD,grad,workingSet);
              
    for (int i=0;i<workingSet.size();i++) { for (int j=0;j<workingSet[i].size();j++) {
    
        step[i][workingSet[i][j]] = 0;
                
        for (int k=0;k<i;k++) { for (int l=0;l<workingSet[k].size();l++) {
                    
            step[i][workingSet[i][j]] -= hess[hindex(k,i,grad.size())][sindex(workingSet[k][l],workingSet[i][j],grad[k].size(),grad[i].size())] * grad[k][workingSet[k][l]];
                        
        } }
        for (int l=0;l<j;l++) {
                    
            step[i][workingSet[i][j]] -= hess[hindex(i,i,grad.size())][hindex(workingSet[i][l],workingSet[i][j],grad[i].size())] * grad[i][workingSet[i][l]];
                        
        }
        for (int l=j;l<workingSet[i].size();l++) {
                    
            step[i][workingSet[i][j]] -= hess[hindex(i,i,grad.size())][hindex(workingSet[i][j],workingSet[i][l],grad[i].size())] * grad[i][workingSet[i][l]];
                        
        }
        for (int k=i+1;k<workingSet.size();k++) { for (int l=0;l<workingSet[k].size();l++) {
                    
            step[i][workingSet[i][j]] -= hess[hindex(i,k,grad.size())][sindex(workingSet[i][j],workingSet[k][l],grad[i].size(),grad[k].size())] * grad[k][workingSet[k][l]];
                        
        } }
                    
    } }

}


// Forward/backward line search for optimal step size

double lineSearch_simple(const Vector &J, const Vector &p, Vector &step, double gamma, const Vector &holdGrad, Vector &grad, Vector &tempJ, Vector &expJ, double S, bool accelerate, double lastAlpha, const IntVector &workingSet, const IntVector &interactingSpins, const IntVector &independentSpins) {
   
    // Try the step with initial alpha
    
    double stepMax     = LInfinity(step, workingSet);
    double alpha       = (stepMax*lastAlpha<MAXSTEPSIZE) ? lastAlpha : MAXSTEPSIZE/stepMax;
    double gradDotStep = innerProduct(grad, step);
    double tempS;
    
    // Sanity check
    if (gradDotStep>0) {
    
        printf("SANITY CHECK FAILED");
        computeDescentStep(holdGrad, sizetolength(J.size()), step, workingSet);
        stepMax     = LInfinity(step, workingSet);
        alpha       = (stepMax<MAXSTEPSIZE) ? 1 : MAXSTEPSIZE/stepMax;
        gradDotStep = innerProduct(holdGrad, step, workingSet);
        
    }
    
    for (int i=0;i<workingSet.size();i++) { for (int j=0;j<workingSet[i].size();j++) {
    
        tempJ[i][workingSet[i][j]] = J[i][workingSet[i][j]] + alpha * step[i][workingSet[i][j]];
        expJ[i][workingSet[i][j]]  = exp(tempJ[i][workingSet[i][j]]);
        
    } }
    
    // If necessary adjust alpha to satisfy weak sufficient decrease
    
    optimizeS_gradOnly(tempJ, expJ, tempS, grad, p, interactingSpins, independentSpins);
    (*regularizeS_gradOnly_ptr)(tempJ, tempS, grad, p, gamma);
    
    int alphaCount = 0;
    bool sufficientDecrease = (tempS<(S + alpha * ARMIJOC * gradDotStep));
    
    while (!sufficientDecrease) {
        
        alphaCount++;
        alpha *= ARMIJOP;
        
        if (stepMax * alpha < MINSTEPSIZE) { alpha = MINSTEPSIZE/stepMax; break; }
            
        for (int i=0;i<workingSet.size();i++) { for (int j=0;j<workingSet[i].size();j++) {
            
            tempJ[i][workingSet[i][j]] = J[i][workingSet[i][j]] + alpha * step[i][workingSet[i][j]];
            expJ[i][workingSet[i][j]]  = exp(tempJ[i][workingSet[i][j]]);
            
        } }
            
        optimizeS_gradOnly(tempJ, expJ, tempS, grad, p, interactingSpins, independentSpins);
        (*regularizeS_gradOnly_ptr)(tempJ, tempS, grad, p, gamma);
        
        if (tempS<(S + alpha * ARMIJOC * gradDotStep)) sufficientDecrease=true;
            
    }
    
    // If sufficient decrease is already satisfied, attempt to step farther
    
    if (alphaCount==0 && accelerate) {
        
        double dotProduct=1;
        double lastS=tempS;
        double holdGradNorm=L2(holdGrad, workingSet);
        
        while (sufficientDecrease) {
            
            alphaCount++;
            
            alpha /= ARMIJOP;
            
            if (stepMax * alpha > MAXSTEPSIZE) { alpha = MAXSTEPSIZE/stepMax; break; }
            
            for (int i=0;i<workingSet.size();i++) { for (int j=0;j<workingSet[i].size();j++) {
                
                tempJ[i][workingSet[i][j]] = J[i][workingSet[i][j]] + alpha * step[i][workingSet[i][j]];
                expJ[i][workingSet[i][j]]  = exp(tempJ[i][workingSet[i][j]]);
                
            } }
            
            optimizeS_gradOnly(tempJ, expJ, tempS, grad, p, interactingSpins, independentSpins);
            (*regularizeS_gradOnly_ptr)(tempJ, tempS, grad, p, gamma);
            
            dotProduct = innerProduct(holdGrad, grad, workingSet) / (holdGradNorm * L2(grad, workingSet));
            
            if ( (tempS>(S + alpha * ARMIJOC * gradDotStep)) || (dotProduct<MAXANGLE) ) { alpha *= ARMIJOP; break; }
            
            lastS=tempS;
            
        }
        
    }
    
    // Return the result
    
    return alpha;
    
}


// Forward/backward line search for optimal step size

double lineSearch_simple(const Vector &J, const Vector &p, Vector &step, double gamma, const Vector &holdGrad, Vector &grad, Vector &tempJ, Vector &expJ, double S, const IntVector &workingSet, const IntVector &interactingSpins, const IntVector &independentSpins) {

    // Try the step with initial alpha
    
    double stepMax     = LInfinity(step, workingSet);
    double alpha       = (stepMax<MAXSTEPSIZE) ? 1 : MAXSTEPSIZE/stepMax;
    double gradDotStep = innerProduct(holdGrad, step, workingSet);
    double tempS;
    
    // Sanity check
    if (gradDotStep>0) {
    
        computeDescentStep(holdGrad, sizetolength(J.size()), step, workingSet);
        stepMax     = LInfinity(step, workingSet);
        alpha       = (stepMax<MAXSTEPSIZE) ? 1 : MAXSTEPSIZE/stepMax;
        gradDotStep = innerProduct(holdGrad, step, workingSet);
        
    }
    
    for (int i=0;i<workingSet.size();i++) { for (int j=0;j<workingSet[i].size();j++) {
    
        tempJ[i][workingSet[i][j]] = J[i][workingSet[i][j]] + alpha * step[i][workingSet[i][j]];
        expJ[i][workingSet[i][j]]  = exp(tempJ[i][workingSet[i][j]]);
        
    } }
    
    // If necessary adjust alpha to satisfy weak sufficient decrease
    
    optimizeS_gradOnly(tempJ, expJ, tempS, grad, p, interactingSpins, independentSpins);
    (*regularizeS_gradOnly_ptr)(tempJ, tempS, grad, p, gamma);
    
    int alphaCount = 0;
    bool sufficientDecrease = (tempS<(S + alpha * ARMIJOC * gradDotStep));
    
    while (!sufficientDecrease) {
        
        alphaCount++;
        alpha *= ARMIJOP;
        
        if (stepMax * alpha < MINSTEPSIZE) { alpha = MINSTEPSIZE/stepMax; break; }
            
        for (int i=0;i<workingSet.size();i++) { for (int j=0;j<workingSet[i].size();j++) {
            
            tempJ[i][workingSet[i][j]] = J[i][workingSet[i][j]] + alpha * step[i][workingSet[i][j]];
            expJ[i][workingSet[i][j]]  = exp(tempJ[i][workingSet[i][j]]);
            
        } }
            
        optimizeS_gradOnly(tempJ, expJ, tempS, grad, p, interactingSpins, independentSpins);
        (*regularizeS_gradOnly_ptr)(tempJ, tempS, grad, p, gamma);
        
        if (tempS<(S + alpha * ARMIJOC * gradDotStep)) sufficientDecrease=true;
            
    }
    
    // If sufficient decrease is already satisfied, attempt to step farther
    
    if (alphaCount==0) {
        
        double dotProduct=1;
        double lastS=tempS;
        double holdGradNorm=L2(holdGrad, workingSet);
        
        while (sufficientDecrease) {
            
            alphaCount++;
            
            alpha /= ARMIJOP;
            
            if (stepMax * alpha > MAXSTEPSIZE) { alpha = MAXSTEPSIZE/stepMax; break; }
            
            for (int i=0;i<workingSet.size();i++) { for (int j=0;j<workingSet[i].size();j++) {
                
                tempJ[i][workingSet[i][j]] = J[i][workingSet[i][j]] + alpha * step[i][workingSet[i][j]];
                expJ[i][workingSet[i][j]]  = exp(tempJ[i][workingSet[i][j]]);
                
            } }
            
            optimizeS_gradOnly(tempJ, expJ, tempS, grad, p, interactingSpins, independentSpins);
            (*regularizeS_gradOnly_ptr)(tempJ, tempS, grad, p, gamma);
            
            dotProduct = innerProduct(holdGrad, grad, workingSet) / (holdGradNorm * L2(grad, workingSet));
            
            if ( (tempS>(S + alpha * ARMIJOC * gradDotStep)) || (dotProduct<MAXANGLE) ) { alpha *= ARMIJOP; break; }
            
            lastS=tempS;
            
        }
        
    }
    
    // Return the result
    
    return alpha;
    
}


// Forward/backward line search for optimal step size

double lineSearch_interp(const Vector &J, const Vector &p, Vector &step, double gamma, const Vector &holdGrad, Vector &grad, Vector &tempJ, Vector &expJ, double S, const IntVector &workingSet, const IntVector &interactingSpins, const IntVector &independentSpins, double armijoS) {

    double c1 = 1.0e-5;
    double c2 = 1.0e+5;
   
    // Try the step with initial alpha
    
    double stepMax     = LInfinity(step, workingSet);
    double alpha       = (stepMax<MAXSTEPSIZE) ? 1 : MAXSTEPSIZE/stepMax;
    double gradDotStep = innerProduct(holdGrad, step, workingSet);
    double stepNorm    = L2(step, workingSet);
    double gradNorm    = L2(grad, workingSet);
    double tempS;
    
    // Sanity check on step direction
    
    if ( (fabs(gradDotStep) < c1 * gradNorm * gradNorm) || (stepNorm > c2 * gradNorm) ) {
    
        computeDescentStep(holdGrad, sizetolength(J.size()), step, workingSet);
        return lineSearch_simple(J, p, step, gamma, holdGrad, grad, tempJ, expJ, S, workingSet, interactingSpins, independentSpins);
        
    }
    
    if (gradDotStep > 0) {
    
        for (int i=0;i<workingSet.size();i++) { for (int j=0;j<workingSet[i].size();j++) step[i][j] = -step[i][j]; }
        
    }
    
    // Test step with maximum alpha
        
    for (int i=0;i<workingSet.size();i++) { for (int j=0;j<workingSet[i].size();j++) {
    
        tempJ[i][workingSet[i][j]] = J[i][workingSet[i][j]] + alpha * step[i][workingSet[i][j]];
        expJ[i][workingSet[i][j]]  = exp(tempJ[i][workingSet[i][j]]);
        
    } }
    
    // If necessary adjust alpha to satisfy weak sufficient decrease
    
    optimizeS_gradOnly(tempJ, expJ, tempS, grad, p, interactingSpins, independentSpins);
    (*regularizeS_gradOnly_ptr)(tempJ, tempS, grad, p, gamma);
    
    int alphaCount = 0;
    bool sufficientDecrease = (tempS<(armijoS + alpha * ARMIJOC * gradDotStep));
    
    while (!sufficientDecrease) {
        
        alphaCount++;
        alpha *= ARMIJOP;
        
        if (stepMax * alpha < MINSTEPSIZE) { alpha = MINSTEPSIZE/stepMax; break; }
            
        for (int i=0;i<workingSet.size();i++) { for (int j=0;j<workingSet[i].size();j++) {
            
            tempJ[i][workingSet[i][j]] = J[i][workingSet[i][j]] + alpha * step[i][workingSet[i][j]];
            expJ[i][workingSet[i][j]]  = exp(tempJ[i][workingSet[i][j]]);
            
        } }
            
        optimizeS_gradOnly(tempJ, expJ, tempS, grad, p, interactingSpins, independentSpins);
        (*regularizeS_gradOnly_ptr)(tempJ, tempS, grad, p, gamma);
        
        if (tempS<(armijoS + alpha * ARMIJOC * gradDotStep)) sufficientDecrease=true;
            
    }
    
    // Return the result
    
    return alpha;
    
}


// Update couplings according to the specified step direction and step size (sparse version)

void makeStep(Vector &J, const Vector &p, const Vector &step, double alpha, double gamma, Vector &grad, Vector &expJ, double &S, const IntVector &workingSet, const IntVector &interactingSpins, const IntVector &independentSpins) {

    for (int i=0;i<workingSet.size();i++) { for (int j=0;j<workingSet[i].size();j++) {
    
        J[i][workingSet[i][j]]    += alpha * step[i][workingSet[i][j]];
        expJ[i][workingSet[i][j]]  = exp(J[i][workingSet[i][j]]);
        
    } }
    
    optimizeS_gradOnly(J, expJ, S, grad, p, interactingSpins, independentSpins);
    (*regularizeS_gradOnly_ptr)(J, S, grad, p, gamma);

}





/*****************************************
       
        M A I N    F U N C T I O N S
        
   Mixed gradient descent and Newton's
   method, L0 + L2-norm regularization,
   see "On learning discrete graphical
   models using greedy methods."
 
*****************************************/





// Compute S and J by minimizing S^*, subject to L0 + L2-norm regularization

void computeSandJ_L0(const Vector &p, int length, double gamma, double gamma0, Vector &J, double &S) {

    // This routine is a controller for the L0 regularization procedure. The idea is to iteratively try minimize the entropy
    // with different sets of couplings according to a backtracking algorithm. gamma0 is the penalty for nonzero couplings.
	
    double nu       = 0.5; 
    int testForward = 3;
    int testBack    = 3;
    
    // workingSet gives the set of variables which are either not regularized or nonzero
    
    IntVector workingSet(J.size(), std::vector<int>());
    
    for (int i=0;i<length;i++)        { for (int j=0;j<J[i].size();j++) {                      workingSet[i].push_back(j); } }
    for (int i=length;i<J.size();i++) { for (int j=0;j<J[i].size();j++) { if (fabs(J[i][j])>0) workingSet[i].push_back(j); } }
    
    // Make temporary variables
    
    Vector expJ(J);
    Vector tempJ(J);
    Vector step(J);
    Vector grad(J);
    Vector holdGrad(J);
    Vector hess;
    std::vector<double> linearHess;
    std::vector<double> inverseC;
    std::vector<double> inverseN;
    std::vector<double> inverseD;
    
    // Get starting entropy
    
    computeSandJ_L2_sparse(p, length, gamma, J, S, workingSet, expJ, tempJ, step, grad, holdGrad, hess, linearHess, inverseC, inverseN, inverseD);
    
    double tryS=S;
    Vector tryJ(J);
    Vector tryGrad(grad);
    
    // Start L0 loop
    
    IntVector addCouplings(testForward,std::vector<int>(2,0));
    IntVector removeCouplings(testBack,std::vector<int>(2,0));
    
    while (true) {
    
        // Get pairs with largest gradients
        
        for (int i=0;i<addCouplings.size();i++) { addCouplings[i][0] = 0; addCouplings[i][1] = 0; }
        
        for (int i=length;i<J.size();i++) { for (int j=0;j<J[i].size();j++) {
        
            for (int k=0;k<addCouplings.size();k++) {
            
                if ( (fabs(grad[i][j])>fabs(grad[addCouplings[k][0]][addCouplings[k][1]]) || addCouplings[k][0]==0) && fabs(J[i][j])==0) {
                
                    for (int l=(int)addCouplings.size()-1;l>k;l--) addCouplings[l] = addCouplings[l-1];
                    
                    addCouplings[k][0]=i;
                    addCouplings[k][1]=j;
                    
                    break;
                    
                }
                
            }
            
        } }
        
        // Add new pairs to workingSet and test couplings
        
        bool acceptCouplings=false;
        
        for (int i=0;i<addCouplings.size();i++) {
        
            // Sanity check
        
            if (addCouplings[i][0]<length)                         { break; }
            if (fabs(J[addCouplings[i][0]][addCouplings[i][1]])>0) { printf("PROBLEM - ATTEMPTED TO ADD A COUPLING THAT EXISTS!\n"); break; }
            
            // Add new coupling and get cluster entropy
        
            insertInPlace(workingSet[addCouplings[i][0]],addCouplings[i][1]);
            computeSandJ_L2_sparse(p, length, gamma, tryJ, tryS, workingSet, expJ, tempJ, step, tryGrad, holdGrad, hess, linearHess, inverseC, inverseN, inverseD);
            
            // If entropy drops significantly, accept step and restart
            
            if (S - tryS > gamma0) {
            
                for (int k=0;k<J.size();k++) { for (int j=0;j<J[k].size();j++) {
                
                    J[k][j]=tryJ[k][j];
                    grad[k][j]=tryGrad[k][j];
                    
                } }
            
                S=tryS;
                
                acceptCouplings=true;
                break;
                
            }
            
            // Else remove the coupling from the set and retry
            
            else {
            
                for (int k=0;k<J.size();k++) { for (int j=0;j<J[k].size();j++) {
                
                    tryJ[k][j]=J[k][j];
                    tryGrad[k][j]=grad[k][j];
                    
                } }
                
                eraseInPlace(workingSet[addCouplings[i][0]],addCouplings[i][1]);
                
            }
        
        }
        
        if (acceptCouplings) continue;
        
        // Get pairs with smallest couplings
        
        for (int i=0;i<removeCouplings.size();i++) { removeCouplings[i][0] = 0; removeCouplings[i][1] = 0; }
        
        for (int i=length;i<workingSet.size();i++) { for (int j=0;j<workingSet[i].size();j++) {
        
            for (int k=0;k<removeCouplings.size();k++) {
            
                if (removeCouplings[k][0]==0) {
                
                    removeCouplings[k][0]=i;
                    //removeCouplings[k][1]=j;
                    removeCouplings[k][1]=workingSet[i][j];
                    break;
                    
                }
            
                if (fabs(J[i][workingSet[i][j]])<fabs(J[removeCouplings[k][0]][removeCouplings[k][1]])) {
                
                    for (int l=(int)removeCouplings.size()-1;l>k;l--) removeCouplings[l] = removeCouplings[l-1];
                    
                    removeCouplings[k][0]=i;
                    //removeCouplings[k][1]=j;
                    removeCouplings[k][1]=workingSet[i][j];
                    break;
                    
                }
                
            }
            
        } }
        
        // Remove pairs from workingSet and test couplings
        
        for (int i=0;i<removeCouplings.size();i++) {
        
            // Sanity check
        
            if (removeCouplings[i][0]<length)                             { break; }
            if (fabs(J[removeCouplings[i][0]][removeCouplings[i][1]])==0) { printf("PROBLEM - ATTEMPTED TO REMOVE A COUPLING THAT IS ZERO!\n"); break; }
            
            // Remove coupling and get cluster entropy
            
            tryJ[removeCouplings[i][0]][removeCouplings[i][1]]=0;
            eraseInPlace(workingSet[removeCouplings[i][0]],removeCouplings[i][1]);
            
            computeSandJ_L2_sparse(p, length, gamma, tryJ, tryS, workingSet, expJ, tempJ, step, tryGrad, holdGrad, hess, linearHess, inverseC, inverseN, inverseD);
            
            // If entropy does not increase significantly, remove coupling and restart
            
            if (tryS - S < nu * gamma0) {
            
                for (int k=0;k<J.size();k++) { for (int j=0;j<J[k].size();j++) {
                
                    J[k][j]=tryJ[k][j];
                    grad[k][j]=tryGrad[k][j];
                    
                } }
                
                S=tryS;
                
                acceptCouplings=true;
                break;
                
            }
            
            // Else return the coupling to workingSet set and retry
            
            else {
            
                insertInPlace(workingSet[removeCouplings[i][0]],removeCouplings[i][1]);
                
                for (int k=0;k<J.size();k++) { for (int j=0;j<J[k].size();j++) {
                
                    tryJ[k][j]=J[k][j];
                    tryGrad[k][j]=grad[k][j];
                    
                } }
                
            }
        
        }
        
        if (acceptCouplings) continue;
        else                 break;
        
    }
    
}


// Compute S and J by minimizing S^*, subject to L2-norm regularization (sparse version)

void computeSandJ_L2_sparse(const Vector &p, int length, double gamma, Vector &J, double &S, const IntVector &workingSet, Vector &expJ, Vector &tempJ, Vector &step, Vector &grad, Vector &holdGrad, Vector &hess, std::vector<double> &linearHess, std::vector<double> &inverseC, std::vector<double> &inverseN, std::vector<double> &inverseD) {
    
    // Define and initialize temporary variables
    
    double eps; // Error
    int    invLength=0;
    bool   descentStep=true;
    double lastAlpha=1.0;
    
    IntVector interactingSpins;
    IntVector independentSpins;
    
    initialize(J, length, expJ, hess, invLength, workingSet, interactingSpins, independentSpins);
    
    linearHess.resize(invLength * (invLength + 1) / 2);
    inverseC.resize(linearHess.size());
    inverseN.resize(linearHess.size());
    inverseD.resize(invLength);
    
    for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) { step[i][j]=0; tempJ[i][j]=0; } }
    
    // Get initial values
    
    optimizeS(J, expJ, S, grad, hess, p, interactingSpins, independentSpins);
    (*regularizeS_ptr)(J, S, grad, hess, p, gamma);
    
    eps = LInfinity(grad, workingSet);
    
    // If not optimal, continue to loop
    
    if (eps>EPSG) descentStep = useDescent(grad, workingSet);
    
    // Create vector of past entropies for relaxed step length choice
    
    std::vector<double> lastS(10,S);
    double armijoS=LInfinity(lastS);
    
    //DEBUG
    int count=0;
    int cNewt=0;
    int cGrad=0;
    
    while (eps>EPSG) {
    
        //DEBUG
//        if (count==0 && length>1) {
//            printf("working set: {\n");
//            for (int i=0;i<workingSet.size();i++) {
//                printf("\t{");
//                for (int j=0;j<workingSet[i].size();j++) printf(" %d",workingSet[i][j]);
//                printf(" }\n");
//            }
//            printf("\t}\n");
//            if (descentStep) printf("Next step: descent\n");
//            else             printf("Next step: Newton\n");
//            printf("J: {\n");
//            for (int i=0;i<J.size();i++) {
//                printf("\t{");
//                for (int j=0;j<J[i].size();j++) printf(" %.4e",J[i][j]);
//                printf(" }\n");
//            }
//            printf("\t}\n");
//            printf("step: {\n");
//            for (int i=0;i<J.size();i++) {
//                printf("\t{");
//                for (int j=0;j<J[i].size();j++) printf(" %.4e",step[i][j]);
//                printf(" }\n");
//            }
//            printf("\t}\n");
//            printf("grad: {\n");
//            for (int i=0;i<J.size();i++) {
//                printf("\t{");
//                for (int j=0;j<J[i].size();j++) printf(" %.4e",grad[i][j]);
//                printf(" }\n");
//            }
//            printf("\t}\n");
//        }
        //DEBUG
        
        //DEBUG
        count++;
        if (descentStep) cGrad++;
        else             cNewt++;
        
        // Compute step direction using gradient descent/Newton's method
        
        if (descentStep) computeDescentStep(grad, length, step, workingSet);
        else             computeNewtonStep(grad, hess, linearHess, inverseC, inverseN, inverseD, step, workingSet);
        
        // Compute optimal step size with weak sufficient decrease condition and make step
        
        for (int i=0;i<workingSet.size();i++) { for (int j=0;j<workingSet[i].size();j++) { holdGrad[i][workingSet[i][j]]=grad[i][workingSet[i][j]]; } }
        
        double alpha;
        if (descentStep) alpha = lineSearch_simple(J, p, step, gamma, holdGrad, grad, tempJ, expJ, armijoS, 1, lastAlpha, workingSet, interactingSpins, independentSpins);
        else             alpha = lineSearch_simple(J, p, step, gamma, holdGrad, grad, tempJ, expJ, armijoS, 0, 1, workingSet, interactingSpins, independentSpins);
        makeStep(J, p, step, alpha, gamma, grad, expJ, S, workingSet, interactingSpins, independentSpins);
        lastAlpha = alpha;
        
        // Update error and choose the next step type
        
        eps         = LInfinity(grad, workingSet);
        descentStep = useDescent(grad, workingSet);
        
        //if (alpha<1e-8 && descentStep) descentStep=false;
        
        if (!descentStep && eps>EPSG) { optimizeS(J, expJ, S, grad, hess, p, interactingSpins, independentSpins); (*regularizeS_ptr)(J, S, grad, hess, p, gamma); }
        
        // Update entropy list for relaxed step length choice
        
        lastS.insert(lastS.begin(),S);
        lastS.pop_back();
        armijoS=LInfinity(lastS);
        
        //DEBUG
//        if (count==10000) {
//            optimizeS(J, expJ, S, grad, hess, p, interactingSpins, independentSpins);
//            (*regularizeS_ptr)(J, S, grad, hess, p, gamma);
//            printf("independentSpins: {\n");
//            for (int i=0;i<independentSpins.size();i++) {
//                printf("\t{");
//                for (int j=0;j<independentSpins[i].size();j++) printf(" %d",independentSpins[i][j]);
//                printf(" }\n");
//            }
//            printf("interactingSpins: {\n");
//            for (int i=0;i<interactingSpins.size();i++) {
//                printf("\t{");
//                for (int j=0;j<interactingSpins[i].size();j++) printf(" %d",interactingSpins[i][j]);
//                printf(" }\n");
//            }
//            printf("p: {\n");
//            for (int i=0;i<p.size();i++) {
//                printf("\t{");
//                for (int j=0;j<p[i].size();j++) printf(" %.8e",p[i][j]);
//                printf(" }\n");
//            }
//            printf("J: {\n");
//            for (int i=0;i<J.size();i++) {
//                printf("\t{");
//                for (int j=0;j<J[i].size();j++) printf(" %.8e",J[i][j]);
//                printf(" }\n");
//            }
//            printf("\t}\n");
//            printf("step: {\n");
//            for (int i=0;i<J.size();i++) {
//                printf("\t{");
//                for (int j=0;j<J[i].size();j++) printf(" %.8e",step[i][j]);
//                printf(" }\n");
//            }
//            printf("\t}\n");
//            printf("grad: {\n");
//            for (int i=0;i<J.size();i++) {
//                printf("\t{");
//                for (int j=0;j<J[i].size();j++) printf(" %.8e",grad[i][j]);
//                printf(" }\n");
//            }
//            printf("\t}\n");
//            optimizeS_gradOnly(J, expJ, S, grad, p, interactingSpins, independentSpins);
//            (*regularizeS_gradOnly_ptr)(J, S, grad, p, gamma);
//            printf("grad (gradOnly): {\n");
//            for (int i=0;i<J.size();i++) {
//                printf("\t{");
//                for (int j=0;j<J[i].size();j++) printf(" %.8e",grad[i][j]);
//                printf(" }\n");
//            }
//            printf("\t}\n");
//        }
//        if ( (count>0 && count<5 && length>1) || (count>10000 && count<10010 && length>1) ) {
//            if (descentStep) printf("Next step: descent, alpha=%.4e\n",alpha);
//            else             printf("Next step: Newton, alpha=%.4e\n",alpha);
//            printf("J: {\n");
//            for (int i=0;i<J.size();i++) {
//                printf("\t{");
//                for (int j=0;j<J[i].size();j++) printf(" %.4e",J[i][j]);
//                printf(" }\n");
//            }
//            printf("\t}\n");
//            printf("step: {\n");
//            for (int i=0;i<J.size();i++) {
//                printf("\t{");
//                for (int j=0;j<J[i].size();j++) printf(" %.4e",step[i][j]);
//                printf(" }\n");
//            }
//            printf("\t}\n");
//            printf("grad: {\n");
//            for (int i=0;i<J.size();i++) {
//                printf("\t{");
//                for (int j=0;j<J[i].size();j++) printf(" %.4e",grad[i][j]);
//                printf(" }\n");
//            }
//            printf("\t}\n");
//        }
        //DEBUG
        
    }
    
    // Return S
    
    optimizeS_gradOnly(J, expJ, S, grad, p, interactingSpins, independentSpins);
    (*regularizeS_ptr)(J, S, grad, hess, p, gamma);
    
}

