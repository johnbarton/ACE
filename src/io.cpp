#include <vector>

#include "io.h"     // Data input and output
#include "tools.h"  // Numerical tools


// Read correlations and number of sites in from a file

void getCorrelations(FILE *input, Vector &p) {

    double c;
    char o;
    
    while (fscanf(input,"%le",&c)==1) {
    
        p.push_back(std::vector<double>());
        (p.back()).push_back(c);
        
        while (fscanf(input,"%c",&o)==1) {
    
            if (o=='\n' || o=='\r') break;
            
            fscanf(input,"%le",&c);
            (p.back()).push_back(c);
            
        }
        
    }
	
}


// Read couplings in from a file

void getCouplings(FILE *input, Vector &p) {

    getCorrelations(input, p);
	
}


// Read secondary structure from file 

void getSS(FILE *input, IntVector &p) {

    int c;
    char o;
    
    while (fscanf(input,"%d",&c)==1) {
    
        p.push_back(std::vector<int>());
        (p.back()).push_back(c);
        
        while (fscanf(input,"%c",&o)==1) {
    
            if (o=='\n' || o=='\r') break;
            
            fscanf(input,"%d",&c);
            (p.back()).push_back(c);
            
        }
    
    }

}


// Read in a multiple sequence alignment in fasta format

void getMSA(FILE *input, std::vector<std::vector<char> > &msa, std::vector<std::string> &tag) {

    char o;
    
    // Find start of sequence recording
    
    while (fscanf(input,"%c",&o)==1) {
    
        if (o=='>') {
        
            msa.push_back(std::vector<char>());
            break;
            
        }
        
    }
    
    while (true) {
            
        // Read in the tag, convert to string, and append
        
        std::vector<char> temptag;
            
        while (fscanf(input,"%c",&o)==1) {
            
            if (o=='\n' || o=='\r') break;
            else temptag.push_back(o);
                
        }
            
        std::string temptag_str(temptag.begin(),temptag.end());
        tag.push_back(temptag_str);
            
        // Read in genetic sequence
        
        bool eof = true;
        
        while (fscanf(input,"%c",&o)==1) {
    
            if (o=='\n' || o=='\r') continue;
            else if (o=='>') {
            
                eof = false;
                msa.push_back(std::vector<char>());
                break;
                
            }
            else (msa.back()).push_back(o);
            
        }
        
        if (eof) break;
    
    }
    
}


// Print MSA to file

void printMSA(FILE *output, const std::vector<std::vector<char> > &msa, const std::vector<std::string> &tag) {

    int linewidth = 100;

    for (int i=0;i<msa.size();i++) {
    
        fprintf(output,">%s",tag[i].c_str());
        
        for (int j=0;j<msa[i].size();j++) {
        
            if (j%linewidth==0) fprintf(output,"\n");
            fprintf(output,"%c",msa[i][j]);
            
        }
        
        fprintf(output,"\n");
        
    }
    
    fflush(output);
    
}


// Print final couplings out to file in a form understood by ising.cpp

void printCouplings(FILE *output, const Vector &J) {
    
    for (int i=0;i<J.size();i++) {
    
        fprintf(output,"%.6e",J[i][0]);
        
        for (int j=1;j<J[i].size();j++) fprintf(output,"\t%.6e",J[i][j]);
        
        fprintf(output,"\n");
        
    }
    
    fflush(output);
	
}


// Print supplementary information to file

void printSupplementaryOutput(FILE *output, double theta, const std::vector<double> &error, double finalS, int maxClusterSize, unsigned long numClusters, unsigned long numSignificantClusters) {
    
    fprintf(output,"%le",theta);
    for (int i=0;i<error.size();i++) fprintf(output,"\t%le",error[i]);
    fprintf(output,"\t%le\t%d\t%lu\t%lu\n",finalS,maxClusterSize,numClusters,numSignificantClusters);
    fflush(output);
	
}


///////////////////////////////


// Read alignment and compute correlations with 3-points

void getAlignment(FILE *input, FILE *weightIn, Vector &J, Vector &p, std::vector<std::vector<std::vector<std::vector<double> > > > &p3, std::vector<double> &pk, std::vector<int> &cons) {
  
    int c = 0;
    int B = 0;
    int N = sizetolength(J.size());
    double Meff = 0;
    std::vector<int> aa;
    std::vector<int> indexi;
    std::vector<int> indexa;
    std::vector<double> weight;

    while (fscanf(input,"%d",&c)==1) { aa.push_back(c); if (c>B) B=c; }
    
    // Get sequence weights from file 
    
    if (weightIn!=NULL) {
        
        getWeights(weightIn,weight);
        for (int i=0;i<weight.size();i++) Meff += weight[i];
        
        printf("Re-weighting: Beff = %.2f with B = %d\n\n",Meff,B);
        
    }
    else {
	
        Meff = (double) B;
        for (int i=0;i<B;i++) weight.push_back(1.0);
        
        printf("No re-weighting vector found: Beff = %.2f with B = %d\n\n",Meff,B);
   
    }
  
    for (int i=0;i<N;i++) { for (int a=0;a<J[i].size();a++) {
      
	  indexi.push_back(i);
	  indexa.push_back(a);
      
    } }

    int cont=-1;
    std::vector<int> mut(B,N);
    std::vector<std::vector<int> > seq(B,std::vector<int>());
    
    for (int i=0;i<aa.size();i++) {
      
        if (aa[i]<0) cont++;
      
        else {
      
            p[indexi[cont]][indexa[cont]] += weight[aa[i]-1]/Meff;
            seq[aa[i]-1].push_back(cont);
            mut[aa[i]-1] -= (int) (indexa[cont]==cons[indexi[cont]]);
        
        }
        
    }
    
    // check for colors outside the alignement
  
    for (int i=0;i<mut.size();i++) { for (int site=0;site<N;site++){
    
        int flag=0;
	  
        for (int j=0;j<seq[i].size();j++) { if (indexi[seq[i][j]]==site) flag=1; }
        
        if (flag==0) mut[i] -= (int) (cons[site]<0);
        
    } }
    
    // compute correlations and p(k) 
    
    for (int i=0;i<B;i++) { for (int j=0;j<seq[i].size();j++) {
    
        int ix = indexi[seq[i][j]];
        int ax = indexa[seq[i][j]];
	  
        for (int k=j+1;k<seq[i].size();k++) {
        
            int jx  = indexi[seq[i][k]];
            int bx  = indexa[seq[i][k]];
            int idx = index(ix,jx,N);
            int sab = sindex(ax,bx,p[ix].size(),p[jx].size());
	      
            p[idx][sab] += weight[i]/Meff;
	      
            for (int l=k+1;l<seq[i].size();l++) {
          
                int kx = indexi[seq[i][l]];
                int cx = indexa[seq[i][l]];
                
                p3[ix][jx][kx][sindex3(ax,bx,cx,p[ix].size(),p[jx].size(),p[kx].size())] += weight[i]/Meff;
	  
            }
        
        }
      
    }
    
    pk[mut[i]] += weight[i]/Meff;
    
    }

}




// Read alignment and compute correlations WITHOUT 3-points

void getAlignment(FILE *input, FILE *weightIn, Vector &J, Vector &p, std::vector<double> &pk, std::vector<int> &cons) {
  
    int c = 0;
    int B = 0;
    int N = sizetolength(J.size());
    double Meff = 0;
    std::vector<int> aa;
    std::vector<int> indexi;
    std::vector<int> indexa;
    std::vector<double> weight;

    while (fscanf(input,"%d",&c)==1) { aa.push_back(c); if (c>B) B=c; }
    
    // Get sequence weights from file 
    
    if (weightIn!=NULL) {
        
        getWeights(weightIn,weight);
        for (int i=0;i<weight.size();i++) Meff += weight[i];
        
        printf("Re-weighting: Beff = %.2f with B = %d\n\n",Meff,B);
        
    }
    else {
	
        Meff = (double) B;
        for (int i=0;i<B;i++) weight.push_back(1.0);
        
        printf("No re-weighting vector found: Beff = %.2f with B = %d\n\n",Meff,B);
   
    }
  
    for (int i=0;i<N;i++) { for (int a=0;a<J[i].size();a++) {
      
	  indexi.push_back(i);
	  indexa.push_back(a);
      
    } }

    int cont=-1;
    std::vector<int> mut(B,N);
    std::vector<std::vector<int> > seq(B,std::vector<int>());
    
    for (int i=0;i<aa.size();i++) {
      
        if (aa[i]<0) cont++;
      
        else {
      
            p[indexi[cont]][indexa[cont]] += weight[aa[i]-1]/Meff;
            seq[aa[i]-1].push_back(cont);
            mut[aa[i]-1] -= (int) (indexa[cont]==cons[indexi[cont]]);
        
        }
        
    }
    
    // check for colors outside the alignement
  
    for (int i=0;i<mut.size();i++) { for (int site=0;site<N;site++){
    
        int flag=0;
	  
        for (int j=0;j<seq[i].size();j++) { if (indexi[seq[i][j]]==site) flag=1; }
        
        if (flag==0) mut[i] -= (int) (cons[site]<0);
        
    } }
    
    // compute correlations and p(k) 
    
    for (int i=0;i<B;i++) { for (int j=0;j<seq[i].size();j++) {
    
        int ix = indexi[seq[i][j]];
        int ax = indexa[seq[i][j]];
	  
        for (int k=j+1;k<seq[i].size();k++) {
        
            int jx  = indexi[seq[i][k]];
            int bx  = indexa[seq[i][k]];
            int idx = index(ix,jx,N);
            int sab = sindex(ax,bx,p[ix].size(),p[jx].size());
	      
            p[idx][sab] += weight[i]/Meff;
	  
        }
      
    }
    
    pk[mut[i]] += weight[i]/Meff;
    
    }

}


// Read consensus from a file

void getConsensus(FILE *input, std::vector<int> &p) {

    int c;
    char o;
    
    while (fscanf(input,"%d",&c)==1) p.push_back(c-1);

}

// Read reweighting vector from a file

void getWeights(FILE *input, std::vector<double> &p) {

    double c;
    char o;
    
    while (fscanf(input,"%le",&c)==1) p.push_back(c);

}

// Print magnetisations

void printMagnetisations(FILE *output, const Vector &J, const Vector &orJ) {
    
    int N = sizetolength(J.size()); 
    
    for (int i=0;i<N;i++) { for (int a=0;a<J[i].size();a++) {
	  
        if (J[i][a]!=0.0 || orJ[i][a]!=0.0) fprintf(output,"%.6e\t%.6e %d %d\n",J[i][a],orJ[i][a],i,a);
        
    } }
    
//    for (int i=0;i<N;i++)  for (int j=0;j<J[i].size();j++) fprintf(output,"%.6e\t%.6e\n",J[i][j],orJ[i][j]); 

    fflush(output);

}


// Print some correlations

void printCorrelations(FILE *outputc, const Vector &Jc, const Vector &orJc, FILE *output, const Vector &J, const Vector &orJ) {
  
    int N    = sizetolength(J.size());
    int cont = 0;
    
      
    for (int i=N;i<J.size();i++) { for (int j=0;j<J[i].size();j++) {
	    
        fprintf(output,"%.6e\t%.6e\n",J[i][j],orJ[i][j]);
	    fprintf(outputc,"%.6e\t%.6e\n",Jc[i][j],orJc[i][j]); 
	
    } }
      
    fflush(output);
	
}

// Print some 3-point correlations

void print3points(FILE *outputc, const std::vector<std::vector<std::vector<std::vector<double> > > > &Jc, const std::vector<std::vector<std::vector<std::vector<double> > > > &orJc,FILE *output, const std::vector<std::vector<std::vector<std::vector<double> > > > &J, const std::vector<std::vector<std::vector<std::vector<double> > > > &orJ,double num) {
   
    int N    = J.size();
    int cont = 0;

    if (num!=0.0) {
      
        printf("Warning: only 3-point correlations bigger than %f will be systematically printed\n", num);
        
        for (int i=0;i<N;i++) { for (int j=i+1;j<N;j++) { for (int k=j+1;k<N;k++) { for (int a=0;a<J[i][j][k].size();a++) {
                
            if (fabs(J[i][j][k][a])>num  || fabs(orJ[i][j][k][a])>num) {
          
                fprintf(output,"%.6e\t%.6e\n",J[i][j][k][a],orJ[i][j][k][a]);
                fprintf(outputc,"%.6e\t%.6e\n",Jc[i][j][k][a],orJc[i][j][k][a]);
            
            }
            else {
                
                if (cont % 50 == 0) {
            
                    fprintf(output,"%.6e\t%.6e\n",J[i][j][k][a],orJ[i][j][k][a]);
                    fprintf(outputc,"%.6e\t%.6e\n",Jc[i][j][k][a],orJc[i][j][k][a]);
            
                }
                
                cont++;
            
            }
            
        } } } }
        
    }
    else {
      
        for (int i=0;i<N;i++) { for (int j=i+1;j<N;j++) { for (int k=j+1;k<N;k++) { for (int l=0;l<J[i][j][k].size();l++) {
          
            fprintf(output,"%.6e\t%.6e\n",J[i][j][k][l],orJ[i][j][k][l]);
            fprintf(outputc,"%.6e\t%.6e\n",Jc[i][j][k][l],orJc[i][j][k][l]);
          
        } } } }
   
    }

    fflush(output);

}


