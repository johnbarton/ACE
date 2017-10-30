#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <vector>
#include <map>
#include <set>

#include "algorithm.h"      // Algorithm declarations
#include "io.h"             // Data input and output
#include "tools.h"          // Numerical tools
#include "inverse.h"        // Computations of cluster entropies
#include "sparseInverse.h"  // Sparse computations of cluster entropies
#include "monteCarlo.h"     // Monte Carlo simulation for checking fit


// GLOBAL VARIABLES

int    N      = 0; // total number of spins
double gamma0 = 0; // strength of the L0 regularization
double gamma2 = 0; // strength of the L2 regularization

int storageSize = 8*sizeof(unsigned long); // number of spins which can be held in a single storage position of the key
int keySize     = 0;                       // number of entries necessary on the key to hold information on all spins

Vector correlations;                                                // the set of all single and pair correlations
std::map<Key,Cluster> *clusterIndex = new std::map<Key,Cluster>();  // the set of all clusters





// Calculate dS and dJ for a given cluster
// This gets the NEGATIVE of dS and dJ without having included S and J (to be input later)
// so that this "guess" for J can be used in the optimization

void computeDSandDJ(int spins[], int clusterSize, Vector &dJ, double &dS) {
    
    // Select the subsets of the cluster and add their contributions to dS and dJ
	
    int subset_mask[clusterSize];
	subset_mask[0]=1;
	for (int i=1;i<clusterSize;i++) subset_mask[i]=0;
    
    Key subsetKey(keySize,0);
    
    unsigned int numSubClusters=pow(2.0,clusterSize)-2;
	
	for (unsigned int n=0;n<numSubClusters;n++) {
		
		// Get the subset (if not found, make a new cluster)
		
        for (int i=0;i<keySize;i++) subsetKey[i]=0;
        
        int subsetSize=0;
        
        for (int i=0;i<clusterSize;i++) { if (subset_mask[i]) {
                
            subsetKey[spins[i]/storageSize] |= (unsigned long) 1 << (spins[i] % storageSize);
            subsetSize++;
                
        } }
        
        //AVOID THIS CHECK IF USING STRICT CONSTRUCTION RULE
        if ((*clusterIndex).count(subsetKey)==0) {
            
            Cluster c(subsetSize);
            makeCluster(c, subsetKey, subsetSize);

            (*clusterIndex)[subsetKey]=c;
            
        }
        
        Cluster const &subCluster=(*clusterIndex)[subsetKey];
		
		// Subtract the contribution to dS

        dS += subCluster.dS;
        
        // Map into the larger cluster and subtract the contribution to dJ
		
		int dJ_count=0;
		
		for (int i=0;i<clusterSize;i++) { if (subset_mask[i]) {
				
            for (int j=0;j<dJ[i].size();j++) dJ[i][j] += subCluster.dJ[dJ_count][j];
            dJ_count++;
				
        } }
        
        for (int i=0;i<clusterSize;i++) { if (subset_mask[i]) {
                
            int off=offset(i,clusterSize);
				
            for (int j=i+1;j<clusterSize;j++) { if (subset_mask[j]) {
						
                for (int k=0;k<dJ[off + j].size();k++) dJ[off + j][k] += subCluster.dJ[dJ_count][k];
                dJ_count++;
						
            } }
				
        } }
        
        // Iterate the subset mask
		
		int i;
		for (i=0;i<clusterSize && subset_mask[i];i++) subset_mask[i]=0;
		if (i<clusterSize) subset_mask[i]=1;
        
	}
    
}


// Make a new cluster of size 1

void makeSingleCluster(Cluster &cluster, const Key &key, int clusterSize) {
    
    // Record individual spins in the key
    
    int spins[clusterSize];
	int n=0;
    
    for (int i=0;i<keySize && n<clusterSize;i++) { for (int j=0;j<storageSize && n<clusterSize;j++) {
            
        if (key[i] & (unsigned long) 1<<j) { spins[n] = storageSize * i + j; n++; }
        
    } }
    
    Vector p(cluster.dJ.size(),std::vector<double>());
    Vector dJ(cluster.dJ.size(),std::vector<double>());
    
	// Get the set of correlations for this cluster and set proper dJ size
	
	for (int i=0;i<clusterSize;i++) {
        
        p[i].resize(correlations[spins[i]].size(),0);
        dJ[i].resize(correlations[spins[i]].size(),0);
        cluster.dJ[i].resize(correlations[spins[i]].size(),0);
        
        for (int j=0;j<correlations[spins[i]].size();j++) p[i][j]=correlations[spins[i]][j];
        
    }
    
	for (int i=0;i<clusterSize-1;i++) { 
        
		int off=offset(spins[i],N);
		
		for (int j=i+1;j<clusterSize;j++) {
			
			p[n].resize(correlations[off + spins[j]].size(),0);
            dJ[n].resize(correlations[off + spins[j]].size(),0);
            cluster.dJ[n].resize(correlations[off + spins[j]].size(),0);
            
            for (int k=0;k<correlations[off + spins[j]].size();k++) p[n][k]=correlations[off + spins[j]][k];
            n++;
			
		}
		
	}
    
    // If the correlations are singular, assign a large field for this site and set dS to zero
    
    double pmax = 0;
    for (int i=0;i<dJ[0].size();i++) { if (p[0][i]>pmax) pmax = p[0][i]; }
    
    if (pmax==0 || pmax==1) {
    
        double hcut = 100;
        for (int i=0;i<cluster.dJ[0].size();i++) {
        
            if (p[0][i]<1) cluster.dJ[0][i] = -hcut;
            else           cluster.dJ[0][i] = 0;
            
        }
        cluster.dS = 0;
        
    }
    
    // Else, compute values and assign them to the cluster
    
    else {
        
        double ptot = 1;
        double dS   = 0;
        for (int i=0;i<dJ[0].size();i++) { dJ[0][i]  = log(p[0][i]); ptot -= p[0][i]; dS -= p[0][i] * log(p[0][i]); }
        for (int i=0;i<dJ[0].size();i++) { dJ[0][i] -= log(ptot);                     dS -= ptot * log(ptot);       }
        
        for (int i=0;i<cluster.dJ.size();i++) { for (int j=0;j<cluster.dJ[i].size();j++) cluster.dJ[i][j]=dJ[i][j]; }
        cluster.dS=dS;
        
    }

}


// Make a new cluster

void makeCluster(Cluster &cluster, const Key &key, int clusterSize) {
    
    // Record individual spins in the key
    
    int spins[clusterSize];
	int n=0;
    
    for (int i=0;i<keySize && n<clusterSize;i++) { for (int j=0;j<storageSize && n<clusterSize;j++) {
            
        if (key[i] & (unsigned long) 1<<j) { spins[n] = storageSize * i + j; n++; }
    
    } }
    
    Vector p(cluster.dJ.size(),std::vector<double>());
    Vector dJ(cluster.dJ.size(),std::vector<double>());
    
	// Get the set of correlations for this cluster and set proper dJ size
	
	for (int i=0;i<clusterSize;i++) {
        
        p[i].resize(correlations[spins[i]].size(),0);
        dJ[i].resize(correlations[spins[i]].size(),0);
        cluster.dJ[i].resize(correlations[spins[i]].size(),0);
        
        for (int j=0;j<correlations[spins[i]].size();j++) p[i][j]=correlations[spins[i]][j];
        
    }
    
	for (int i=0;i<clusterSize-1;i++) { 
        
		int off=offset(spins[i],N);
		
		for (int j=i+1;j<clusterSize;j++) {
			
			p[n].resize(correlations[off + spins[j]].size(),0);
            dJ[n].resize(correlations[off + spins[j]].size(),0);
            cluster.dJ[n].resize(correlations[off + spins[j]].size(),0);
            
            for (int k=0;k<correlations[off + spins[j]].size();k++) p[n][k]=correlations[off + spins[j]][k];
            n++;
			
		}
		
	}
    
//    //DEBUG
//    if (clusterSize>0) {
//    printf("cluster {%d",spins[0]);
//    for (int i=1;i<clusterSize;i++) printf(", %d",spins[i]);
//    printf("}\n"); }
//    //DEBUG
    
	// Compute values and assign them to the cluster
    
    double dS = 0;
    
    (*computeS0andJ0_ptr)(p,clusterSize,gamma2,dJ,dS);
    if (clusterSize>1) computeDSandDJ(spins,clusterSize,dJ,dS);
    
    Vector J(dJ);
    double S=0;
    
    (*computeSandJ_ptr)(p,clusterSize,gamma2,gamma0,J,S);
    for (int i=0;i<cluster.dJ.size();i++) { for (int j=0;j<cluster.dJ[i].size();j++) cluster.dJ[i][j]=(J[i][j]-dJ[i][j]); }
    cluster.dS=S-dS;

}


// Run the cluster selection and creation algorithm

void selectClusters(std::set<Key> &clusters, int clusterSize, int spinCutSize, int pairCutSize, int cutSize, double theta, bool useVerbose, FILE *output) {
	
	std::vector<Key> significantClusters;
    significantClusters.reserve(clusters.size());
    
    Significant isSignificant(theta);
    
    for (std::set<Key>::iterator i=clusters.begin();i!=clusters.end();++i) {
        
        if ( isSignificant((*clusterIndex)[*i]) ) {
	  
            significantClusters.push_back(*i);
            
            // Print cluster for further analysis
            
            if (output!=NULL) {
	  
                int spins[clusterSize];
                int n=0;
                Key key=*i;
    
                for (int i=0;i<keySize && n<clusterSize;i++) { for (int j=0;j<storageSize && n<clusterSize;j++) {
            
                    if (key[i] & (unsigned long) 1<<j) { spins[n] = storageSize * i + j; n++; }
    
                } }
	  
                fprintf(output," %d",spins[0]);
                for (int i=1;i<clusterSize;i++) fprintf(output," %d",spins[i]);
                fprintf(output,"\n");
            
            }
	
        }
        
	}
    
    clusters.clear();
    
    if (useVerbose) printf("Found %d significant cluster(s). ",(int) significantClusters.size());
    
    
    // Build list of potential new cluster elements
    
    std::map<int,int> spinMap;
    std::set<int> spinSet;
    
    std::map<std::vector<int>,int> pairMap;
    std::set<std::vector<int> > pairSet;
    std::vector<int> pair(2,0);
    
    for (int i=0;i<significantClusters.size();i++) {
    
        // Extract list of spins
    
        int spins[clusterSize];
        int n=0;
        
        for (int j=0;j<keySize && n<clusterSize;j++) { for (int k=0;k<storageSize && n<clusterSize;k++) {
            
            if (significantClusters[i][j] & (unsigned long) 1<<k) { spins[n] = storageSize * j + k; n++; }
        
        } }
        
        // Count singles and pairs
        
        for (int j=0;j<clusterSize;j++) {
        
            if (spinMap.count(spins[j])==0) spinMap[spins[j]]  = 1;
            else                            spinMap[spins[j]] += 1;
        
            pair[0]=spins[j];
            
            for (int k=j+1;k<clusterSize;k++) {
            
                pair[1]=spins[k];
                if (pairMap.count(pair)==0) pairMap[pair]  = 1;
                else                        pairMap[pair] += 1;
                
            }
            
        }
    
    }
    
    // Select potential singles
    
    for (std::map<int,int>::iterator i=spinMap.begin();i!=spinMap.end();++i) {
    
           if ((*i).second>=spinCutSize) spinSet.insert(spinSet.end(),(*i).first);
        
    }
    if (useVerbose) printf("%ld spins. ",spinSet.size());
    spinMap.clear();
    
    // Select potential pairs
    
    for (std::map<std::vector<int>,int>::iterator i=pairMap.begin();i!=pairMap.end();++i) {
    
           if ((*i).second>=pairCutSize) pairSet.insert(pairSet.end(),(*i).first);
        
    }
    if (useVerbose) printf("%ld pairs. ",pairSet.size());
    pairMap.clear();

    
    // Find supersets of significant clusters and reassign them to clusters

    Key supersetKey(keySize,0);
    
    std::map<Key,int> buildNum;
    
    for (int i=0;i<significantClusters.size();i++) {
    
        // Extract list of spins, skip if single or pair is absent
    
        int spins[clusterSize];
        int n=0;
        bool halt=false;
    
        for (int j=0;j<keySize && n<clusterSize;j++) { for (int k=0;k<storageSize && n<clusterSize;k++) {
            
            if (significantClusters[i][j] & (unsigned long) 1<<k) { spins[n] = storageSize * j + k; n++; }
            
        } }
        
        for (int j=0;j<clusterSize && !halt;j++) {
        
            if (spinSet.count(spins[j])==0) { halt=true; break; }
        
            pair[0]=spins[j];
            
            for (int k=j+1;k<clusterSize && !halt;k++) {
            
                pair[1]=spins[k];
                if (pairSet.count(pair)==0) { halt=true; break; }
                
            }
            
        }
        
        if (halt) continue;
        
        // Add potential new clusters
        
        for (int j=0;j<keySize;j++) supersetKey[j]=significantClusters[i][j];
		
        for (std::set<int>::iterator j=spinSet.begin();j!=spinSet.end();++j) {
            
            // If cluster is unchanged, skip
            
            if (supersetKey[(*j)/storageSize] & (unsigned long) 1 << ((*j) % storageSize)) continue;
            
            // Test for potential cluster
            
            bool isPossible=true;
            
            for (int k=0;k<clusterSize;k++) {
            
                if (spins[k]<(*j)) { pair[0]=spins[k]; pair[1]=(*j);     }
                else               { pair[0]=(*j);     pair[1]=spins[k]; }
                
                if (pairSet.count(pair)==0) { isPossible=false; break; }
            
            }
            
            // Add to count
            
            if (isPossible || clusterSize==2) {
            
                supersetKey[(*j)/storageSize] |= (unsigned long) 1 << ((*j) % storageSize);

                if (buildNum.count(supersetKey)==0) buildNum[supersetKey]=1;
                else                                buildNum[supersetKey]+=1;
            
                supersetKey[(*j)/storageSize]  = significantClusters[i][(*j)/storageSize];
                
            }
        
		}
		
	}
    
    
    // Cut clusters of of size>3 if not enough subclusters are significant
    
    for (std::map<Key,int>::iterator i=buildNum.begin();i!=buildNum.end();++i) {
    
           if ((*i).second>cutSize) clusters.insert(clusters.end(),(*i).first);
        
    }

}


// Build the list of clusters

void getClusters(std::set<Key> &clusters, int &maxClusterSize, int kmax, double theta, bool lax, bool useVerbose, const std::vector<int> &allowedSites,FILE *output) {

    // Form all allowed pair clusters
    
    if (allowedSites.size()==0) {

        for (int i=0;i<N;i++) { for (int j=i+1;j<N;j++) {
                
            Key clusterKey(keySize,0);
                
            clusterKey[i/storageSize] |= (unsigned long) 1 << (i % storageSize);
            clusterKey[j/storageSize] |= (unsigned long) 1 << (j % storageSize);
                
            clusters.insert(clusterKey);
            
        } }
        
    }
    
    else {
    
        for (int i=0;i<allowedSites.size();i++) { for (int j=i+1;j<allowedSites.size();j++) {
                
            Key clusterKey(keySize,0);
                
            clusterKey[allowedSites[i]/storageSize] |= (unsigned long) 1 << (allowedSites[i] % storageSize);
            clusterKey[allowedSites[j]/storageSize] |= (unsigned long) 1 << (allowedSites[j] % storageSize);
                
            clusters.insert(clusterKey);
            
        } }
        
    }
    
    // Continue to general algorithm
	
	while (clusters.size()>0 && maxClusterSize!=kmax) {
        
        maxClusterSize++;
        if (useVerbose) printf("Computing all significant clusters of size %d. ",maxClusterSize);

        for (std::set<Key>::iterator i=clusters.begin();i!=clusters.end();++i) {
        
            if ((*clusterIndex).count(*i)==0) {
            
                Cluster c(maxClusterSize);
                makeCluster(c,*i,maxClusterSize);
                (*clusterIndex)[*i]=c;
                
            }
            
        }
        
        // Set cutoffs for fast cluster construction
        
        int spinCutSize = maxClusterSize;
        int pairCutSize = maxClusterSize-1;
        int cutSize     = maxClusterSize;
        
        if (lax) {
        
            spinCutSize = 1;
            pairCutSize = 0;
            cutSize     = 2;
            
        }
        
        if (maxClusterSize==2) { spinCutSize = 1; pairCutSize = 0; cutSize = 2; }
        
        selectClusters(clusters,maxClusterSize,spinCutSize,pairCutSize,cutSize,theta,useVerbose,output);
        
        if (useVerbose) printf("Formed %d new cluster(s).\n",(int) clusters.size());
        
    }

}


// Iterate through the map of clusters to obtain final couplings

void getCouplings(Vector &finalJ, double &finalS, unsigned long &numSignificantClusters, double theta) {

    Key subsetKey(keySize,0);

    for (std::map<Key,Cluster>::iterator i=(*clusterIndex).begin();i!=(*clusterIndex).end();++i) {
        
        Significant isSignificant(theta);
		
		if ( isSignificant((*i).second) ) {
        
            // Map cluster to the whole system
            
            int clusterSize = sizetolength((*i).second.dJ.size());
            int spins[clusterSize];
            
            int n=0;
            
            for (int j=0;j<keySize && n<clusterSize;j++) { for (int k=0;k<storageSize && n<clusterSize;k++) {
                    
                if ((*i).first[j] & (unsigned long) 1<<k) { spins[n] = storageSize * j + k; n++; }
                
            } }
        
            // Check if cluster was not previously selected
            
            if (!(*i).second.selected) {
            
                (*i).second.selected = true;
        
                // Iterate through subclusters
                
                unsigned int numSubClusters = pow(2.0,clusterSize) - 2;
                
                int subset_mask[clusterSize];
                subset_mask[0]=1;
                for (int j=1;j<clusterSize;j++) subset_mask[j] = 0;
	
                for (unsigned int n=0;n<numSubClusters;n++) {
		
                    for (int j=0;j<keySize;j++) subsetKey[j] = 0;
                
                    int subsetSize = 0;
                    for (int j=0;j<clusterSize;j++) { if (subset_mask[j]) {
                
                        subsetKey[spins[j]/storageSize] |= (unsigned long) 1 << (spins[j] % storageSize);
                        subsetSize++;
                
                    } }
                    
                    // Mark as superSelected
                    
                    (*clusterIndex)[subsetKey].superSelected = true;
                    
                    // Iterate the subset mask
		
                    int j;
                    for (j=0;j<clusterSize && subset_mask[j];j++) subset_mask[j] = 0;
                    if (j<clusterSize) subset_mask[j]=1;
                    
                }
                
            }
            
            // Count the number of significant clusters
            
            numSignificantClusters++;
            
            // Add the cluster's contribution to dS
            
            finalS += (*i).second.dS;
            
            // Add the contribution to dJ
            
            for (int j=0;j<clusterSize;j++) {
            
                for (int k=0;k<finalJ[spins[j]].size();k++) finalJ[spins[j]][k] += (*i).second.dJ[j][k];
                
            }
            
            for (int j=0;j<clusterSize-1;j++) {
                
                int off=offset(spins[j], N);
                
                for (int k=j+1;k<clusterSize;k++) {
                
                    for (int l=0;l<finalJ[off + spins[k]].size();l++) finalJ[off + spins[k]][l] += (*i).second.dJ[n][l];
                    n++;
                    
                }
				
			}
			
		}
        
	}

}


// Runs the main program

int run(RunParameters &r) {
    
    // Retrieve correlations from file and set system, key sizes
    
    if (FILE *datain = fopen(r.getCorrelationsInfile().c_str(),"r")) {
        
        getCorrelations(datain, correlations);
        
        N = sizetolength(correlations.size());
        keySize = (N + storageSize - 1) / storageSize;
    
        fclose(datain);
    
    }
    else {
    
        printf("Problem retrieving data from file %s! The file may not exist, or it may be inaccessible\n",r.getCorrelationsInfile().c_str());
        return EXIT_FAILURE;
        
    }
    
    // Open supplementary output file
    
    FILE *supout = fopen(r.getSupplementaryOutfile().c_str(),"w");
    
    // If number of data points given, compute optimal regularization strength
    
    if (r.sampleB>0 && r.computeGamma) {
        
        r.gamma0=computeGamma_L0(correlations, r.sampleB);
        r.gamma2=computeGamma_L2(correlations, r.sampleB);
        
    }
    
    // Set kmax to the system size, if no specific maximum size is given 
    
    if (r.kmax<=0) r.kmax=N;
    
    // Sanity checks on theta and gamma input
    
    double theta = r.thetaMax;
    
    if (r.thetaMin<=0) { RunParameters rDefault; r.thetaMin  = rDefault.thetaMin;  }
    if (r.thetaMin>theta)      theta      = r.thetaMin;
    if (r.thetaMax<r.thetaMin) r.thetaMax = r.thetaMin;
    if (r.thetaStep<1) { RunParameters rDefault; r.thetaStep = rDefault.thetaStep; }
    if (r.gamma0<0)    { RunParameters rDefault; r.gamma0    = rDefault.gamma0;    }
    if (r.gamma2<0)    { RunParameters rDefault; r.gamma2    = rDefault.gamma2;    }
    
    gamma0=r.gamma0;
    gamma2=r.gamma2;
    
    bool isBinary = true;
    for (int i=0;i<correlations.size();i++) {
        
        if (correlations[i].size()>1) { isBinary = false; break; }
        
    }
    
    if (r.useRef && isBinary) computeS0andJ0_ptr=&computeS0andJ0_L2_binary;
    else if (r.useRef)        computeS0andJ0_ptr=&computeS0andJ0_L2;
    else                      computeS0andJ0_ptr=&computeS0andJ0_Empty;
    
    if (r.useSparse) computeSandJ_ptr=&computeSandJ_L0;
    else             computeSandJ_ptr=&computeSandJ_L2;
    
    if (r.useGI) {
    
        regularizeS_ptr          = &regularizeS_L2_GI;
        regularizeS_gradOnly_ptr = &regularizeS_L2_gradOnly_GI;
        epsilonP2_ptr            = &epsilonP2_GI;
        epsilonC_ptr             = &epsilonC_GI;
        getMaxError_ptr          = &getMaxError_GI;
        
    }
    else {
    
        regularizeS_ptr          = &regularizeS_L2;
        regularizeS_gradOnly_ptr = &regularizeS_L2_gradOnly;
        epsilonP2_ptr            = &epsilonP2;
        epsilonC_ptr             = &epsilonC;
        getMaxError_ptr          = &getMaxError;
        
    }
    
	if (r.useVerbose) {
    
        printf("Inferring Ising model couplings using theta = %.8e, \t gamma0 = %.8e, \t gamma2 = %.8e...\n",theta,gamma0,gamma2);
        if (r.useSparse) printf("L0 regularization is enabled\n");
        printf("Minimum cluster size = %d, maximum cluster size = %d\n",r.kmin,r.kmax);
        printf("Storage size = %d, key size = %d\n\n",storageSize,keySize);
        
    }
    
    // Read default starting clusters

    std::vector<int> ssSites;

    if (r.useCmap || r.inputClusters) {

        if (r.useVerbose) printf("Selecting clusters from the input list\n");

        IntVector ss;

        if (FILE *ssin = fopen(r.getSecStructInfile().c_str(),"r")) {
            
            getSS(ssin, ss);
            
            int maxsize=0;

            for (int i=0;i<ss.size();i++) {

                Key clusterKey(keySize,0);

                for (int j=0;j<ss[i].size();j++) {

                    clusterKey[ss[i][j]/storageSize] |= (unsigned long) 1 << (ss[i][j] % storageSize);
                    insertInPlace(ssSites,ss[i][j]);

                }

                Cluster c((int)ss[i].size());
                makeCluster(c,clusterKey,(int) ss[i].size());
                c.selected = true;
                (*clusterIndex)[clusterKey]=c;

            }
            
            fclose(ssin);
            
    	}
        
        else printf("Problem retrieving data from file %s! The file may not exist, or it may be inaccessible\n",r.getSecStructInfile().c_str());
        
	}
    
    // MAIN ALGORITHM
    
    // Fix the first record value for theta
    double ThetaRec=theta/r.recordStep;
    std::string filename;
    
    // Clusters of size one are treated specially
	
	for (int i=0;i<N;i++) {
		
        Key clusterKey(keySize,0);
        
        clusterKey[i/storageSize] |= (unsigned long) 1 << (i % storageSize);
        
		Cluster c(1);
        makeCluster(c,clusterKey,1);
        (*clusterIndex)[clusterKey]=c;
        
	}
    
    if (r.useVerbose) printf("Computing all clusters of size 1. Found %d cluster(s).\n", N);
    
    // Construct clusters only using sites with non-singular correlations (if ssSites is empty)
    
    if (ssSites.size()==0) { for (int i=0;i<N;i++) {
    
        Key clusterKey(keySize,0);
    
        clusterKey[i/storageSize] |= (unsigned long) 1 << (i % storageSize);
        
        if ((*clusterIndex)[clusterKey].dS!=0) ssSites.push_back(i);
        
    } }
    
    int maxClusterSize=1;
    std::set<Key> clusters;
    FILE *output = NULL;
    if (r.recClusters) output = fopen(r.getClusterOutfile().c_str(),"w");
    if (r.inputClusters==false) getClusters(clusters, maxClusterSize, r.kmax, theta, r.useLax, r.useVerbose, ssSites, output);
    if (r.recClusters) fclose(output);
	
	// GET FINAL COUPLINGS
    
    if (r.useVerbose) printf("Found all significant clusters at threshold %.8e. Getting final couplings.\n",theta);
	
    Vector finalJ0(correlations.size(),std::vector<double>());
    for (int i=0;i<correlations.size();i++) finalJ0[i].resize(correlations[i].size(),0);
    double finalS0=0;
    (*computeS0andJ0_ptr)(correlations,N,gamma2,finalJ0,finalS0);
    
    Vector finalJ(finalJ0);
    double finalS=finalS0;
    
    unsigned long numSignificantClusters=0;
    unsigned long numClusters=(*clusterIndex).size();
    
    getCouplings(finalJ, finalS, numSignificantClusters, theta);
    
    // CHECK ERROR
    
    std::vector<double> error(3,0);
    
    unsigned long lastNumSignificantClusters=numSignificantClusters;
    
    getError(correlations, finalJ, N, r.sampleB, r.b, r.runs, gamma2, ALPHA, error);
    
    // While errors are > 1, lower threshold and loop
    
    // Initialize variables for theta recording
    Vector trecJ(finalJ);
    std::vector<double> minError(3,1000);    
    double recTheta=theta;
    
    while ( (error[0]>1 || error[1]>1 || error[2]>1 || maxClusterSize<r.kmin) && maxClusterSize<r.kmax && theta>r.thetaMin && r.inputClusters==false) {
        
        // First record data
        
        FILE *infout = fopen(r.getCouplingsOutfile().c_str(),"w");
        printCouplings(infout, finalJ);
        fclose(infout);
        
        printSupplementaryOutput(supout, theta, error, finalS, maxClusterSize, numClusters, numSignificantClusters);
        
        if (r.recClusterCover) {
        
            FILE *cout = fopen(r.getClusterCoverOutfile().c_str(),"w");

            for (std::map<Key,Cluster>::iterator i=(*clusterIndex).begin();i!=(*clusterIndex).end();++i) { if ((*i).second.selected && !(*i).second.superSelected) {

                // Map cluster to the whole system

                int clusterSize = sizetolength((*i).second.dJ.size());
                int spins[clusterSize];

                int n = 0;

                for (int j=0;j<keySize && n<clusterSize;j++) { for (int k=0;k<storageSize && n<clusterSize;k++) {

                    if ((*i).first[j] & (unsigned long) 1<<k) { spins[n] = storageSize * j + k; n++; }

                } }

                // Print cluster content

                fprintf(cout,"%.6e\t",(*i).second.dS);
                for (int j=0;j<clusterSize;j++) fprintf(cout," %d",spins[j]);
                fprintf(cout,"\n");
                fflush(cout);

            } }
            
            fclose(cout);
            
        }
        
        // Find min error on correlations and record J and h
        
        if (error[2]<minError[2] && maxClusterSize>=r.kmin) {
        
            minError[0] = error[0];
            minError[1] = error[1];
            minError[2] = error[2];
        
            for (int i=0;i<finalJ.size();i++) { for (int j=0;j<finalJ[i].size();j++) trecJ[i][j]=finalJ[i][j]; }
            recTheta=theta;
            
        }
	
        // Write down J and h with min epsilonMax if theta<ThetaRec
	
        if (theta<ThetaRec && r.recordStep!=0) {
	  
            // Update ThetaRec
            ThetaRec=theta/r.recordStep;
	    
            char value[64];
            sprintf(value, "%.6f",recTheta);
	    
            std::string filename = r.getOutfile_TH()+"_th"+value+".j";
            FILE *infout = fopen(filename.c_str(),"w");
            printCouplings(infout,trecJ);
            fclose(infout);
	    
            filename = r.getOutfile_TH()+".th";
            FILE *out = fopen(filename.c_str(),"a");
            fprintf(out,"%.6f  %d\n",recTheta,1000);
            fclose(out);
	    
            minError[0] = 1000;
            minError[1] = 1000;
            minError[2] = 1000;
            
        }
	
        // Lower threshold and rerun steps above but do not compute cluster if it already exists
        
        theta /= r.thetaStep;
        maxClusterSize=1;
        clusters.clear();
        
        FILE *output = NULL;
        if (r.recClusters) output = fopen(r.getClusterOutfile().c_str(),"w");
        getClusters(clusters, maxClusterSize, r.kmax, theta, r.useLax, r.useVerbose, ssSites, output);
        if (r.recClusters) fclose(output);
        
        // Get final couplings
        
        if (r.useVerbose) printf("Found all significant clusters at threshold %.8e. Getting final couplings.\n",theta);
        
        finalS=finalS0;
        for (int i=0;i<finalJ.size();i++) { for (int j=0;j<finalJ[i].size();j++) finalJ[i][j]=finalJ0[i][j]; }
        
        numSignificantClusters=0;
        numClusters=(*clusterIndex).size();
        
        getCouplings(finalJ, finalS, numSignificantClusters, theta);
        
        // Compute error again if the number of significant clusters has changed
        
        if (numSignificantClusters!=lastNumSignificantClusters) getError(correlations, finalJ, N, r.sampleB, r.b, r.runs, gamma2, ALPHA, error);
        
        lastNumSignificantClusters=numSignificantClusters;
        
    }
    
    // Record final data and exit
    
    FILE *infout = fopen(r.getCouplingsOutfile().c_str(),"w");
    printCouplings(infout,finalJ);
    fclose(infout);

    printSupplementaryOutput(supout, theta, error, finalS, maxClusterSize, numClusters, numSignificantClusters);
    fclose(supout);
    
    return EXIT_SUCCESS;

}

