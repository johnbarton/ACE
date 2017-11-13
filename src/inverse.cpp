#include <vector>

#include "tools.h"      // Numerical tools
#include "inverse.h"    // Computations of cluster entropies


//#define EPSG 1.0e-10
//#define MAXANGLE 0.998        
//#define MAXSTEPSIZE 5.0
//#define MINSTEPSIZE 1.0e-10
//#define NEWTONGRADSIZE  3.0e-3
//#define NEWTONGRADRATIO 1.0e-8
//#define ARMIJOC 1.0e-2
//#define ARMIJOP 0.5
//#define STEMCUT 1.0e-14

#define EPSG 1.0e-9
#define MAXANGLE 0.998        
#define MAXSTEPSIZE 5.0
#define MINSTEPSIZE 1.0e-11
#define NEWTONGRADSIZE  3.0e-3
#define NEWTONGRADRATIO 1.0e-8
#define ARMIJOC 1.0e-5
#define ARMIJOP 0.5
#define STEMCUT 1.0e-14

// GLOBAL VARIABLES

void (*computeS0andJ0_ptr)(const Vector &, int, double, Vector &, double &);
void (*computeSandJ_ptr)(const Vector &, int, double, double, double, Vector &, double &);
void (*regularizeS_ptr)(const Vector &, double &, Vector &, Vector &, const Vector &, double, double);
void (*regularizeS_gradOnly_ptr)(const Vector &, double &, Vector &, const Vector &, double, double);





/*****************************************
       
   P A R T I T I O N    F U N C T I O N
 
*****************************************/





// A reference entropy call that simply returns 0 for all couplings, fields, and the entropy

void computeS0andJ0_Empty(const Vector &p, int length, double unused, Vector &J0, double &S0) {
    
    for (int i=0;i<J0.size();i++) { for (int j=0;j<J0[i].size();j++) J0[i][j]=0; }
    S0=0;
    
}


// Recursive computation of the partition function

double getZ(int L, const Vector &expJ, std::vector<int> &bt, std::vector<int> &spins, Vector &grad, Vector &hess, int length, double stem) {
    
    double temp = 1.0;
        
    // Compute correlations
    
    for (int j=(int)bt.size()-1;j>=0;j--) {
    
        grad[bt[j]][spins[j]] += stem;
        hess[hindex(bt[j],bt[j],expJ.size())][hindex(spins[j],spins[j],expJ[bt[j]].size())] += stem;
            
        int offj = offset(bt[j],length);
            
        for (int k=j-1;k>=0;k--) {
            
            int sjk = sindex(spins[j],spins[k],expJ[bt[j]].size(),expJ[bt[k]].size());
            
            grad[offj + bt[k]][sjk] += stem;
                
            hess[hindex(bt[j],bt[k],expJ.size())][sjk] += stem;
            hess[hindex(bt[j],offj+bt[k],expJ.size())][sindex(spins[j],sjk,expJ[bt[j]].size(),expJ[offj+bt[k]].size())] += stem;
            hess[hindex(bt[k],offj+bt[k],expJ.size())][sindex(spins[k],sjk,expJ[bt[k]].size(),expJ[offj+bt[k]].size())] += stem;
            hess[hindex(offj+bt[k],offj+bt[k],expJ.size())][hindex(sjk,sjk,expJ[offj+bt[k]].size())] += stem;

            int offk = offset(bt[k],length);
                
            for (int l=k-1;l>=0;l--) {
                
                int sjl = sindex(spins[j],spins[l],expJ[bt[j]].size(),expJ[bt[l]].size());
                int skl = sindex(spins[k],spins[l],expJ[bt[k]].size(),expJ[bt[l]].size());
                
                hess[hindex(bt[j],offk+bt[l],expJ.size())][sindex(spins[j],skl,expJ[bt[j]].size(),expJ[offk+bt[l]].size())] += stem;
                hess[hindex(bt[k],offj+bt[l],expJ.size())][sindex(spins[k],sjl,expJ[bt[k]].size(),expJ[offj+bt[l]].size())] += stem;
                hess[hindex(bt[l],offj+bt[k],expJ.size())][sindex(spins[l],sjk,expJ[bt[l]].size(),expJ[offj+bt[k]].size())] += stem;
                    
                hess[hindex(offj+bt[k],offj+bt[l],expJ.size())][sindex(sjk,sjl,expJ[offj+bt[k]].size(),expJ[offj+bt[l]].size())] += stem;
                hess[hindex(offj+bt[k],offk+bt[l],expJ.size())][sindex(sjk,skl,expJ[offj+bt[k]].size(),expJ[offk+bt[l]].size())] += stem;
                hess[hindex(offj+bt[l],offk+bt[l],expJ.size())][sindex(sjl,skl,expJ[offj+bt[l]].size(),expJ[offk+bt[l]].size())] += stem;
                    
                int offl=offset(bt[l],length);
                    
                for (int m=l-1;m>=0;m--) {
                    
                    int sjm = sindex(spins[j],spins[m],expJ[bt[j]].size(),expJ[bt[m]].size());
                    int skm = sindex(spins[k],spins[m],expJ[bt[k]].size(),expJ[bt[m]].size());
                    int slm = sindex(spins[l],spins[m],expJ[bt[l]].size(),expJ[bt[m]].size());
                        
                    hess[hindex(offj+bt[k],offl+bt[m],expJ.size())][sindex(sjk,slm,expJ[offj+bt[k]].size(),expJ[offl+bt[m]].size())] += stem;
                    hess[hindex(offj+bt[l],offk+bt[m],expJ.size())][sindex(sjl,skm,expJ[offj+bt[l]].size(),expJ[offk+bt[m]].size())] += stem;
                    hess[hindex(offj+bt[m],offk+bt[l],expJ.size())][sindex(sjm,skl,expJ[offj+bt[m]].size(),expJ[offk+bt[l]].size())] += stem;
                    
                }
                
            }
            
        }
        
    }
    
    // Truncation condition
    
    //if (stem<STEMCUT) return 1.0;
    
    // Branch further down the tree
    
    for(int i = L-1; i >= 0; i--) {
        
        bt.push_back(i);
        
        int off=offset(i,length);
        
        for (int p=0;p<expJ[i].size();p++) {
            
            double prodJ = expJ[i][p];
            
            for (int j=0;j<bt.size()-1;j++) prodJ *= expJ[off+bt[j]][sindex(p,spins[j],expJ[i].size(),expJ[bt[j]].size())];
            
            spins.push_back(p);
            
            temp += prodJ * getZ(i, expJ, bt, spins, grad, hess, length, prodJ * stem);
            
            spins.pop_back();
        
        }
        
        bt.pop_back();
    
    }
    
    // Pass up the tree
    
    return temp;

}


// Recursive computation of the cross entropy

void optimizeS(const Vector &J, const Vector &expJ, double &S, Vector &grad, Vector &hess, const Vector &p) {

    // Set variables
    
    int length=(int) sizetolength(J.size());

    for (int i=0;i<grad.size();i++) { for (int j=0;j<grad[i].size();j++) grad[i][j]=0; }
    for (int i=0;i<hess.size();i++) { for (int j=0;j<hess[i].size();j++) hess[i][j]=0; }
    
    // Get partition function and correlations
    
    std::vector<int> bt;
    std::vector<int> spins;

    double Z = getZ(length, expJ, bt, spins, grad, hess, length, 1.0);
    
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

double getZ_gradOnly(int L, const Vector &expJ, std::vector<int> &bt, std::vector<int> &spins, Vector &grad, int length, double stem) {
    
    double temp = 1.0;
        
    // Compute correlations
    
    for (int j=(int)bt.size()-1;j>=0;j--) {
    
        grad[bt[j]][spins[j]] += stem;
            
        int offj = offset(bt[j],length);
            
        for (int k=j-1;k>=0;k--) grad[offj + bt[k]][sindex(spins[j],spins[k],expJ[bt[j]].size(),expJ[bt[k]].size())] += stem;
        
    }
    
    // Truncation condition
    
    //if (stem<STEMCUT) return 1.0;
    
    // Branch further down the tree
    
    for(int i = L-1; i >= 0; i--) {
        
        bt.push_back(i);
        
        int off=offset(i,length);
        
        for (int p=0;p<expJ[i].size();p++) {
            
            double prodJ = expJ[i][p];
            
            for (int j=0;j<bt.size()-1;j++) prodJ *= expJ[off+bt[j]][sindex(p,spins[j],expJ[i].size(),expJ[bt[j]].size())];
            
            spins.push_back(p);
            
            temp += prodJ * getZ_gradOnly(i, expJ, bt, spins, grad, length, prodJ * stem);
            
            spins.pop_back();
        
        }
        
        bt.pop_back();
    
    }
    
    // Pass up the tree
    
    return temp;

}


// Recursive computation of the cross entropy (gradient only)

void optimizeS_gradOnly(const Vector &J, const Vector &expJ, double &S, Vector &grad, const Vector &p) {

    // Set variables
    
    int length=(int) sizetolength(J.size());

    for (int i=0;i<grad.size();i++) { for (int j=0;j<grad[i].size();j++) grad[i][j]=0; }
    
    // Get partition function and correlations
    
    std::vector<int> bt;
    std::vector<int> spins;

    double Z = getZ_gradOnly(length, expJ, bt, spins, grad, length, 1.0);
    
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

void initialize(const Vector &J, int length, Vector &expJ, Vector &hess, int &invLength) {

    invLength=0;

    for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) expJ[i][j] = exp(J[i][j]); }
    
    hess.resize( (J.size() * (J.size() + 1)) / 2 );
    
    for (int i=0;i<length;i++) {
    
        hess[hindex(i,i,J.size())].resize( (J[i].size() * (J[i].size() + 1)) / 2, 0);                      // Resize field-field diagonal
        
        for (int j=i+1;j<length;j++)      hess[hindex(i,j,J.size())].resize(J[i].size() * J[j].size(), 0); // Resize field-field off-diagonal
        for (int j=length;j<J.size();j++) hess[hindex(i,j,J.size())].resize(J[i].size() * J[j].size(), 0); // Resize field-coupling
        
        invLength += (int) J[i].size();
        
    }
    for (int i=length;i<J.size();i++) {
    
        hess[hindex(i,i,J.size())].resize( (J[i].size() * (J[i].size() + 1)) / 2, 0);                      // Resize coupling-coupling diagonal
        
        for (int j=i+1;j<J.size();j++)    hess[hindex(i,j,J.size())].resize(J[i].size() * J[j].size(), 0); // Resize coupling-coupling off-diagonal
        
        invLength += (int) J[i].size();
        
    }

}


// Choose to take a step using Newton method or gradient descent

bool useDescent(const Vector &grad) {
    
    double gradMax=0;
    double gradMin=1;
        
    for (int i=0;i<grad.size();i++) { for (int j=0;j<grad[i].size();j++) {
            
        if (fabs(grad[i][j])>gradMax)                       gradMax=fabs(grad[i][j]);
        if (fabs(grad[i][j])<gradMin && fabs(grad[i][j])>0) gradMin=fabs(grad[i][j]);
            
    } }
        
    if (gradMin<EPSG) gradMin=EPSG;
    
    return ( gradMax>NEWTONGRADSIZE || (gradMin/gradMax)<NEWTONGRADRATIO );

}


// Compute step size using gradient descent

void computeDescentStep(const Vector &grad, int length, Vector &step) {

    double totalGrad=0;
    if (length==1) totalGrad=1.0;
        
    for (int i=0;i<grad.size();i++) { for (int j=0;j<grad[i].size();j++) { totalGrad += fabs(grad[i][j]); } }
    for (int i=0;i<grad.size();i++) { for (int j=0;j<grad[i].size();j++) { step[i][j] = -grad[i][j]/totalGrad; } }

}


// Compute step size using Newton's method

void computeNewtonStep(const Vector &grad, Vector &hess, std::vector<double> &linearHess, int length, int invLength, Vector &step) {

    linearHess.resize(invLength * (invLength + 1) / 2);
    modifiedCholeskyInverse(hess,linearHess,invLength,grad);
    
    for (int i=0;i<grad.size();i++) { for (int j=0;j<grad[i].size();j++) {
    
        step[i][j] = 0;
                
        for (int k=0;k<i;k++) {             for (int l=0;l<grad[k].size();l++) step[i][j] -= hess[hindex(k,i,grad.size())][sindex(l,j,grad[k].size(),grad[i].size())] * grad[k][l]; }
        for (int l=0;l<j;l++) {                                                step[i][j] -= hess[hindex(i,i,grad.size())][hindex(l,j,grad[i].size())]                * grad[i][l]; }
        for (int l=j;l<grad[i].size();l++) {                                   step[i][j] -= hess[hindex(i,i,grad.size())][hindex(j,l,grad[i].size())]                * grad[i][l]; }
        for (int k=i+1;k<grad.size();k++) { for (int l=0;l<grad[k].size();l++) step[i][j] -= hess[hindex(i,k,grad.size())][sindex(j,l,grad[i].size(),grad[k].size())] * grad[k][l]; }
        
    } }

}


// Forward/backward line search for optimal step size

double lineSearch_simple(const Vector &J, const Vector &p, Vector &step, double gamma, double gammah, const Vector &holdGrad, Vector &grad, Vector &tempJ, Vector &expJ, double S, bool accelerate, double lastAlpha) {
   
    // Try the step with initial alpha
    
    double stepMax     = LInfinity(step);
    double alpha       = (stepMax*lastAlpha<MAXSTEPSIZE) ? lastAlpha : MAXSTEPSIZE/stepMax;
    double gradDotStep = innerProduct(grad, step);
    double tempS;
        
    for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) {
    
        tempJ[i][j] = J[i][j] + alpha * step[i][j];
        expJ[i][j]  = exp(tempJ[i][j]);
        
    } }
    
    // If necessary adjust alpha to satisfy weak sufficient decrease
    
    optimizeS_gradOnly(tempJ, expJ, tempS, grad, p);
    (*regularizeS_gradOnly_ptr)(tempJ, tempS, grad, p, gamma, gammah);
    
    int alphaCount = 0;
    bool sufficientDecrease = (tempS<(S + alpha * ARMIJOC * gradDotStep));
    
    while (!sufficientDecrease) {
        
        alphaCount++;
        alpha *= ARMIJOP;
            
        if (stepMax * alpha < MINSTEPSIZE) { alpha = MINSTEPSIZE/stepMax; break; }
            
        for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) {
        
            tempJ[i][j] = J[i][j] + alpha * step[i][j];
            expJ[i][j] = exp(tempJ[i][j]);
            
        } }
            
        optimizeS_gradOnly(tempJ, expJ, tempS, grad, p);
        (*regularizeS_gradOnly_ptr)(tempJ, tempS, grad, p, gamma, gammah);
        
        if (tempS<(S + alpha * ARMIJOC * gradDotStep)) sufficientDecrease=true;
            
    }
    
    // If sufficient decrease is already satisfied, attempt to step farther
    
    if (alphaCount==0 && accelerate) {
        
        double dotProduct=1;
        double lastS=tempS;
        double holdGradNorm=L2(holdGrad);
        
        while (sufficientDecrease) {
            
            alphaCount++;
            
            alpha /= ARMIJOP;
            
            if (stepMax * alpha > MAXSTEPSIZE) { alpha = MAXSTEPSIZE/stepMax; break; }
            
            for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) {
            
                tempJ[i][j] = J[i][j] + alpha * step[i][j];
                expJ[i][j] = exp(tempJ[i][j]);
                
            } }
            
            optimizeS_gradOnly(tempJ, expJ, tempS, grad, p);
            (*regularizeS_gradOnly_ptr)(tempJ, tempS, grad, p, gamma, gammah);
            
            dotProduct = innerProduct(holdGrad, grad) / (holdGradNorm * L2(grad));
            
            if ( (tempS>(S + alpha * ARMIJOC * gradDotStep)) || (dotProduct<MAXANGLE) ) { alpha *= ARMIJOP; break; }
            
            lastS=tempS;
            
        }
        
    }
    
    // Return the result
    
    return alpha;
    
}


// Forward/backward line search for optimal step size

double lineSearch_interp(const Vector &J, const Vector &p, Vector &step, double gamma, double gammah, const Vector &holdGrad, Vector &grad, Vector &tempJ, Vector &expJ, double S) {

    double c1 = 1.0e-5;
    double c2 = 1.0e+5;
   
    // Try the step with initial alpha
    
    double stepMax     = LInfinity(step);
    double alpha       = (stepMax<MAXSTEPSIZE) ? 1 : MAXSTEPSIZE/stepMax;
    double gradDotStep = innerProduct(grad, step);
    double stepNorm    = L2(step);
    double gradNorm    = L2(grad);
    double tempS;
    
    // Sanity check on step direction
    
    if ( (fabs(gradDotStep) < c1 * gradNorm * gradNorm) || (stepNorm > c2 * gradNorm) ) {// || (gradDotStep > 0) ) {
    
        computeDescentStep(holdGrad, sizetolength(J.size()), step);
        return lineSearch_simple(J, p, step, gamma, gammah, holdGrad, grad, tempJ, expJ, S, 0, alpha);
        
    }
    
    if (gradDotStep > 0) {
    
        for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) step[i][j] = -step[i][j]; }
        
    }
    
    // Try the step with initial alpha
    
    for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) {
    
        tempJ[i][j] = J[i][j] + alpha * step[i][j];
        expJ[i][j]  = exp(tempJ[i][j]);
        
    } }
    
    // If necessary adjust alpha to satisfy weak sufficient decrease
    
    optimizeS_gradOnly(tempJ, expJ, tempS, grad, p);
    (*regularizeS_gradOnly_ptr)(tempJ, tempS, grad, p, gamma, gammah);
    
    int alphaCount = 0;
    bool sufficientDecrease = (tempS<(S + alpha * ARMIJOC * gradDotStep));
    
    if (!sufficientDecrease) {
    
        double lastAlpha=alpha;
        double lastdSdJ=gradDotStep/stepNorm;
        double lastS=S;
        
        // Update alpha and make trial step
        
        alpha = -lastdSdJ * alpha * alpha / ( 2 * (tempS - lastS - alpha * lastdSdJ) );
        
        for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) {
        
            tempJ[i][j] = J[i][j] + alpha * step[i][j];
            expJ[i][j] = exp(tempJ[i][j]);
            
        } }
        
        optimizeS_gradOnly(tempJ, expJ, tempS, grad, p);
        (*regularizeS_gradOnly_ptr)(tempJ, tempS, grad, p, gamma, gammah);
        
        if (tempS<(S + alpha * ARMIJOC * gradDotStep)) sufficientDecrease=true;
        
        // Loop if sufficient decrease fails
        
        while (!sufficientDecrease && alpha>MAXSTEPSIZE/stepMax && alpha<1.0) {
        
            alphaCount++;
            
            // Compute step variables
            
            double dSdJ=innerProduct(grad,step)/stepNorm;
            double d1=dSdJ+lastdSdJ-3*(lastS-tempS)/(lastAlpha-alpha);
            double d2=sign(alpha-lastAlpha)*sqrt(d1*d1-dSdJ*lastdSdJ);
            
            // Update alpha and make trial step
            
            double lastAlphaHold=alpha;
            
            alpha = alpha - (alpha - lastAlpha) * (dSdJ + d2 - d1) / (dSdJ - lastdSdJ + 2 * d2);
            
            lastdSdJ=dSdJ;
            lastAlpha=lastAlphaHold;
            lastS=tempS;
            
            for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) {
                
                tempJ[i][j] = J[i][j] + alpha * step[i][j];
                expJ[i][j] = exp(tempJ[i][j]);
                
            } }
            
            optimizeS_gradOnly(tempJ, expJ, tempS, grad, p);
            (*regularizeS_gradOnly_ptr)(tempJ, tempS, grad, p, gamma, gammah);
            
            if (tempS<(S + alpha * ARMIJOC * gradDotStep)) sufficientDecrease=true;
        
        }
        
    }
    
    // If sufficient decrease is already satisfied, attempt to step farther (accelerate)
    
    if (alphaCount==0) {
        
        double dotProduct=1;
        double lastS=tempS;
        double holdGradNorm=L2(holdGrad);
        
        while (sufficientDecrease) {
            
            alphaCount++;
            
            alpha /= ARMIJOP;
            
            if (stepMax * alpha > MAXSTEPSIZE) { alpha = MAXSTEPSIZE/stepMax; break; }
            
            for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) {
            
                tempJ[i][j] = J[i][j] + alpha * step[i][j];
                expJ[i][j] = exp(tempJ[i][j]);
                
            } }
            
            optimizeS_gradOnly(tempJ, expJ, tempS, grad, p);
            (*regularizeS_gradOnly_ptr)(tempJ, tempS, grad, p, gamma, gammah);
            
            dotProduct = innerProduct(holdGrad, grad) / (holdGradNorm * L2(grad));
            
            if ( (tempS>(S + alpha * ARMIJOC * gradDotStep)) || (dotProduct<MAXANGLE) ) { alpha *= ARMIJOP; break; }
            
            lastS=tempS;
            
        }
        
    }
    
    // Return the result
    
    return alpha;
    
}



// Update couplings according to the specified step direction and step size

void makeStep(Vector &J, const Vector &p, const Vector &step, double alpha, double gamma, double gammah, Vector &grad, Vector &expJ, double &S) {

    for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) {
    
        J[i][j]    += alpha * step[i][j];
        expJ[i][j]  = exp(J[i][j]);
        
    } }
            
    optimizeS_gradOnly(J, expJ, S, grad, p);
    (*regularizeS_gradOnly_ptr)(J, S, grad, p, gamma, gammah);

}





/*****************************************
       
       M A I N    F U N C T I O N S
        
   Mixed gradient descent and Newton's
   method, L2-norm regularization.
 
*****************************************/





// Compute S0 and J0 given a set of correlations with L2 regularization, BINARY VARIABLES ONLY

void computeS0andJ0_L2_binary(const Vector &p, int length, double gamma, Vector &J0, double &S0) {

    // Make the covariance matrix
    
    double c[(length*(length+1))/2];
	
	for (int i=0;i<length;i++) {
		
		c[hindex(i,i,length)] = 1;
        
        int off = offset(i,length);
		
		for (int j=i+1;j<length;j++) {
            
            if (p[i][0]>0 && p[j][0]>0) c[hindex(i,j,length)] = (p[off + j][0] - p[i][0] * p[j][0]) / sqrt(p[i][0] * (1 - p[i][0]) * p[j][0] * (1 - p[j][0]));
            else                        c[hindex(i,j,length)] = 0;
            
        }

		
	}
    
    // Get eigenvalues and eigenvectors
    
    std::vector<std::vector<double> > eigenvectors(length);
    double eigenvalues[length];
    
    for (int i=0;i<length;i++) { eigenvectors[i].resize(length,0); eigenvectors[i][i]=1; }
    
    if (length==2) {
        
        eigenvalues[0] = 1+c[1];
        eigenvalues[1] = 1-c[1];

        eigenvectors[0][0] = 1/sqrt(2);
        eigenvectors[0][1] = 1/sqrt(2);
        eigenvectors[1][0] = 1/sqrt(2);
        eigenvectors[1][1] = -1/sqrt(2);
        
    }
    
    else symmetricQR(c, length, eigenvalues, eigenvectors);
    
    // Compute S0
    
    double mhat[length];
    S0 = 0;
    
    for (int i=0;i<length;i++) {
        
        mhat[i] = ((eigenvalues[i] - gamma) + sqrt(pow((eigenvalues[i] - gamma), 2.0) + 4 * gamma)) / 2.0;
        S0     += (log(mhat[i]) + 1 - mhat[i]);
        
    }
    
    S0 /= 2;
    
    // Compute the J' = Q^T Lambda Q matrix
    
    for (int i=0;i<length;i++) {
        
        int off = offset(i,length);
        
        for (int j=i+1;j<length;j++) {
            
            J0[off + j][0] = 0;
            
            for (int k=0;k<length;k++) J0[off + j][0] += eigenvectors[i][k] * eigenvectors[j][k] * (1 - 1 / mhat[k]);
            
            J0[off + j][0] /= sqrt(p[i][0] * (1 - p[i][0]) * p[j][0] * (1 - p[j][0]));
            
        }
        
    }
    
    // Compute h0 from standard approach, using J0 already calculated
    
    for (int i=0;i<length;i++) {
        
        double temp = 0;
        
        for        (int j=0;j<i;j++) temp += J0[index(j,i,length)][0] * ( (p[index(j,i,length)][0] - p[i][0] * p[j][0]) * (p[i][0] - 0.5) / (p[i][0] * (1 - p[i][0])) - p[j][0] );
        for (int j=i+1;j<length;j++) temp += J0[index(i,j,length)][0] * ( (p[index(i,j,length)][0] - p[i][0] * p[j][0]) * (p[i][0] - 0.5) / (p[i][0] * (1 - p[i][0])) - p[j][0] );
        
        J0[i][0] = temp;
        
    }
    
}


// Compute S0 and J0 given a set of correlations with L2 regularization

void computeS0andJ0_L2(const Vector &p, int length, double gamma, Vector &J0, double &S0) {

    // Make auxiliary variables
    
    int Neff = 0;
    for (int i=0;i<length;i++) Neff+=(int) p[i].size();

    // Make the covariance matrix
    
    std::vector<double> c((Neff*(Neff+1))/2, 0);
	int count = 0;
    
	for (int i=0;i<length;i++) { for (int a=0;a<p[i].size();a++) {
		
		c[hindex(count,count,Neff)] = 1;
        
        int count2 = count+1;
        
        // Sum over states at the same site
        
        for (int b=a+1;b<p[i].size();b++) {
        
            if (p[i][a]>0 && p[i][b]>0) c[hindex(count,count2,Neff)] = -(p[i][a] * p[i][b]) / sqrt(p[i][a] * (1 - p[i][a]) * p[i][b] * (1 - p[i][b]));
            else                        c[hindex(count,count2,Neff)] = 0;
            
            count2++;
            
        }
		
        // Sum over states at other sites
        
		for (int j=i+1;j<length;j++) { for (int b=0;b<p[j].size();b++) {
        
            int idx  = index(i,j,length);
            int sidx = sindex(a,b,p[i].size(),p[j].size());
            
            if (p[i][a]>0 && p[j][b]>0) c[hindex(count,count2,Neff)] = (p[idx][sidx] - (p[i][a] * p[j][b])) / sqrt(p[i][a] * (1 - p[i][a]) * p[j][b] * (1 - p[j][b]));
            else                        c[hindex(count,count2,Neff)] = 0;
            
            count2++;
            
        } }

        count++;
		
	} }
    
    // Get eigenvalues and eigenvectors
    
    std::vector<std::vector<double> > eigenvectors(Neff);
    double eigenvalues[Neff];
    
    for (int i=0;i<Neff;i++) { eigenvectors[i].resize(Neff,0); eigenvectors[i][i]=1; }
    
    if (Neff==2) {
        
        eigenvalues[0] = 1+c[1];
        eigenvalues[1] = 1-c[1];

        eigenvectors[0][0] = 1/sqrt(2);
        eigenvectors[0][1] = 1/sqrt(2);
        eigenvectors[1][0] = 1/sqrt(2);
        eigenvectors[1][1] = -1/sqrt(2);
        
    }
    
    else symmetricQR(c, Neff, eigenvalues, eigenvectors);
    
    // Compute S0
    
    std::vector<double> mhat(Neff,0);
    S0 = 0;
    
    for (int i=0;i<Neff;i++) {
        
        mhat[i] = ((eigenvalues[i] - gamma) + sqrt(pow((eigenvalues[i] - gamma), 2.0) + 4 * gamma)) / 2.0;
        S0     += (log(mhat[i]) + 1 - mhat[i]);
        
    }
    
    S0 /= 2;
    
    // Compute the J' = Q^T Lambda Q matrix
    
    std::vector<std::vector<double> > Jdiag(length);
    count = 0;
    
    for (int i=0;i<length;i++) Jdiag[i].resize((p[i].size()*(p[i].size()-1))/2,0);

    for (int i=0;i<length;i++) { for (int a=0;a<p[i].size();a++) {
    
        int count2 = count+1;
        
        // Sum over states at the same site (store diagonal Js)
        
        for (int b=a+1;b<p[i].size();b++) {
        
            int idx = index(a,b,p[i].size())-p[i].size();

            Jdiag[i][idx] = 0;

            for (int k=0;k<Neff;k++) Jdiag[i][idx] += eigenvectors[count][k] * eigenvectors[count2][k] * (1 - 1 / mhat[k]);
            
            Jdiag[i][idx] /= sqrt(p[i][a] * (1 - p[i][a]) * p[i][b] * (1 - p[i][b]));
            
            count2++;
            
        }
        
        // Sum over states at other sites
        
        for (int j=i+1;j<length;j++) { for (int b=0;b<p[j].size();b++) {
            
            int idx  = index(i,j,length);
            int sidx = sindex(a,b,p[i].size(),p[j].size());
            
            J0[idx][sidx] = 0;
            
            for (int k=0;k<Neff;k++) J0[idx][sidx] += eigenvectors[count][k] * eigenvectors[count2][k] * (1 - 1 / mhat[k]);
            
            J0[idx][sidx] /= sqrt(p[i][a] * (1 - p[i][a]) * p[j][b] * (1 - p[j][b]));
            
            count2++;
            
        } }
        
        count++;
        
    } }
    
    // Compute h0 from standard approach, using J0 already calculated
    
    for (int i=0;i<length;i++) { for (int a=0;a<p[i].size();a++) {
        
        double temp = 0;
        
        // diagonal
        for (int b=0;b<a;b++)             temp += Jdiag[i][index(b,a,p[i].size())-p[i].size()] * ( -(p[i][a] * p[i][b]) * (p[i][a] - 0.5) / (p[i][a] * (1 - p[i][a])) - p[i][b] );
        for (int b=a+1;b<p[i].size();b++) temp += Jdiag[i][index(a,b,p[i].size())-p[i].size()] * ( -(p[i][a] * p[i][b]) * (p[i][a] - 0.5) / (p[i][a] * (1 - p[i][a])) - p[i][b] );
        
        // off-diagonal
        for (int j=0;j<i;j++) { for (int b=0;b<p[j].size();b++) {
        
            int idx  = index(j,i,length);
            int sidx = sindex(b,a,p[j].size(),p[i].size());
        
            temp += J0[idx][sidx] * ( (p[idx][sidx] - (p[i][a] * p[j][b])) * (p[i][a] - 0.5) / (p[i][a] * (1 - p[i][a])) - p[j][b] );
            
        } }
        for (int j=i+1;j<length;j++) { for (int b=0;b<p[j].size();b++) {
        
            int idx  = index(i,j,length);
            int sidx = sindex(a,b,p[i].size(),p[j].size());
        
            temp += J0[idx][sidx] * ( (p[idx][sidx] - (p[i][a] * p[j][b])) * (p[i][a] - 0.5) / (p[i][a] * (1 - p[i][a])) - p[j][b] );
        
        } }
        
        J0[i][a] = temp;
        
    } }
    
}


// Compute S and J by minimizing S^*, subject to L2-norm regularization

void computeSandJ_L2(const Vector &p, int length, double gamma, double gammah, double unused, Vector &J, double &S) {

    // Define and initialize temporary variables
    
    double eps;                 // Error
    int    invLength=0;
    bool   descentStep=true;
    double lastAlpha=1.0;
    Vector expJ(J);
    Vector tempJ(J);
    Vector step(J);
    Vector tempStep(J);
    Vector grad(J);
    Vector hess;
    std::vector<double> linearHess;
    
    initialize(J, length, expJ, hess, invLength);
    
    // Do first iteration manually to set J
    
    optimizeS(J, expJ, S, grad, hess, p);
    (*regularizeS_ptr)(J, S, grad, hess, p, gamma, gammah);
    
    Vector holdGrad(grad);
    eps = LInfinity(grad);
    
    // Create vector of past entropies for relaxed step length choice
    
    std::vector<double> lastS(10,S);
    double armijoS=LInfinity(lastS);
    
//    ///////
//    //DEBUG
//    if (length>0) {
//        printf("J: {\n");
//        for (int i=0;i<J.size();i++) {
//            printf("\t{");
//            for (int j=0;j<J[i].size();j++) printf(" %.4e",J[i][j]);
//            printf(" }\n");
//        }
//        printf("\t}\n");
//        printf("grad: {\n");
//        for (int i=0;i<J.size();i++) {
//            printf("\t{");
//            for (int j=0;j<J[i].size();j++) printf(" %.4e",grad[i][j]);
//            printf(" }\n");
//        }
//        printf("\t}\n");
////        printf("hess: {\n");
////        for (int i=0;i<hess.size();i++) {
////            printf("\t{");
////            for (int j=0;j<hess[i].size();j++) printf(" %.4e",hess[i][j]);
////            printf(" }\n");
////        }
////        printf("\t}\n");
//    }
//    //DEBUG
//    ///////
    
    // If not optimal, continue to loop
    
    if (eps>EPSG) descentStep = useDescent(grad);
    
    int count = 0;
    while (eps>EPSG) {
        
        if (descentStep) computeDescentStep(grad, length, step);
        else             computeNewtonStep(grad, hess, linearHess, length, invLength, step);
        
        // Compute optimal step size with weak sufficient decrease condition and make step
        
        for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) holdGrad[i][j]=grad[i][j]; }
        
        double alpha;
        if (descentStep) alpha = lineSearch_simple(J, p, step, gamma, gammah, holdGrad, grad, tempJ, expJ, armijoS, 1, lastAlpha);
        else             alpha = lineSearch_simple(J, p, step, gamma, gammah, holdGrad, grad, tempJ, expJ, armijoS, 0, 1);
        makeStep(J, p, step, alpha, gamma, gammah, grad, expJ, S);
        lastAlpha = alpha;
        
        // Update error and choose the next step type
        
        eps         = LInfinity(grad);
        descentStep = useDescent(grad);
        
        if (!descentStep && eps>EPSG) { optimizeS(J, expJ, S, grad, hess, p); (*regularizeS_ptr)(J, S, grad, hess, p, gamma, gammah); }
        
        // Update entropy list for relaxed step length choice
        // Switch out of relaxed step choice if optimization takes a long time
        
        if (count>1000) armijoS = S;
        else {
            lastS.insert(lastS.begin(),S);
            lastS.pop_back();
            armijoS=LInfinity(lastS);
        }
        count++;
        
//        ///////
//        //DEBUG
//        if (count>1000 && count<1010) {
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
////            printf("hess: {\n");
////            for (int i=0;i<hess.size();i++) {
////                printf("\t{");
////                for (int j=0;j<hess[i].size();j++) printf(" %.4e",hess[i][j]);
////                printf(" }\n");
////            }
////            printf("\t}\n");
//        }
//        //DEBUG
//        ///////
        
        
    }
    
    //DEBUG
    //printf("iterations: %d\n",count);
    //DEBUG
    
    // Return S
    
    optimizeS_gradOnly(J, expJ, S, grad, p);
    (*regularizeS_ptr)(J, S, grad, hess, p, gamma, gammah);
    
}


// Modifies the supplied model entropy given a set of couplings, with L2 regularization

void regularizeS_L2(const Vector &J, double &S, Vector &grad, Vector &hess, const Vector &p, double gamma, double gammah) {
    
    int length=(int) sizetolength(J.size());
    double gamma_J_squared=0;
    
    // Regularization for fields
    
    for (int i=0;i<length;i++) { for (int j=0;j<J[i].size();j++) {
            
        gamma_J_squared += gammah * J[i][j] * J[i][j];
        grad[i][j]      += 2 * gammah * gamma * J[i][j];
    
        hess[hindex(i,i,J.size())][hindex(j,j,J[i].size())] += 2 * gammah * gamma;
            
    } }
    
    // Regularization for couplings
    
    for (int i=length;i<J.size();i++) { for (int j=0;j<J[i].size();j++) {
            
        gamma_J_squared += J[i][j] * J[i][j];
        grad[i][j]      += 2 * gamma * J[i][j];
    
        hess[hindex(i,i,J.size())][hindex(j,j,J[i].size())] += 2 * gamma;
            
    } }
	
	// Return the result
	
    gamma_J_squared *= gamma;
	S               += gamma_J_squared;
	
}


// Modifies the supplied model entropy given a set of couplings, with L2 regularization

void regularizeS_L2_gradOnly(const Vector &J, double &S, Vector &grad, const Vector &p, double gamma, double gammah) {
    
    int length=(int) sizetolength(J.size());
    double gamma_J_squared=0;
    
    // Regularization for fields
    
    for (int i=0;i<length;i++) { for (int j=0;j<J[i].size();j++) {
            
        gamma_J_squared += gammah * J[i][j] * J[i][j];
        grad[i][j]      += 2 * gammah * gamma * J[i][j];
    
    } }

    // Regularization for couplings
    
    for (int i=length;i<J.size();i++) { for (int j=0;j<J[i].size();j++) {
            
        gamma_J_squared += J[i][j] * J[i][j];
        grad[i][j]      += 2 * gamma * J[i][j];
    
    } }

	// Return the result
    
    gamma_J_squared *= gamma;
	S               += gamma_J_squared;
	
}

// Modifies the supplied model entropy given a set of couplings, with L2 regularization

void regularizeS_L2_GI(const Vector &J, double &S, Vector &grad, Vector &hess, const Vector &p, double gamma, double gammah) {
    
    int length=(int) sizetolength(J.size());
    double gamma_J_squared=0;
    
    // Regularization for fields
    
    for (int i=0;i<length;i++) { for (int j=0;j<J[i].size();j++) {
            
        gamma_J_squared += gammah * J[i][j] * J[i][j];
        grad[i][j]      += 2 * gammah * gamma * J[i][j];
    
        hess[hindex(i,i,J.size())][hindex(j,j,J[i].size())] += 2 * gammah * gamma;
            
    } }

    // GAUGE INVARIANT Regularization for couplings
    
    for (int i=0;i<length;i++) { for (int j=i+1;j<length;j++) {
    
        // Get quantities for computing K_ij^ab
    
        int Jidx = index(i,j,length);
        int hidx = hindex(Jidx,Jidx,J.size());
        double inv_qi = 1/((double)J[i].size()+1);
        double inv_qj = 1/((double)J[j].size()+1);
    
        double suma[J[i].size()];
        double sumb[J[j].size()];
        double sumall=0;
        
        for (int a=0;a<J[i].size();a++) suma[a]=0;
        for (int b=0;b<J[j].size();b++) sumb[b]=0;
        
        for (int a=0;a<J[i].size();a++) { for (int b=0;b<J[j].size();b++) {
            
            suma[a] += J[Jidx][sindex(a,b,J[i].size(),J[j].size())];
            sumb[b] += J[Jidx][sindex(a,b,J[i].size(),J[j].size())];
            
        } }
        
        for (int a=0;a<J[i].size();a++) sumall += suma[a];
        
        // Compute K_ij^ab and contribution to regularization term and Hessian
        
        for (int a=0;a<J[i].size();a++) { for (int b=0;b<J[j].size();b++) {
        
            int sab  = sindex(a,b,J[i].size(),J[j].size());
            double K = J[Jidx][sab] - (suma[a] * inv_qj) - (sumb[b] * inv_qi) + (sumall * inv_qi * inv_qj);
        
            gamma_J_squared                            += K * K;
            grad[Jidx][sab]                            += 2 * gamma * K;
            hess[hidx][hindex(sab,sab,J[Jidx].size())] += 2 * gamma;
            
            for (int c=a;c<J[i].size();c++) hess[hidx][hindex(sab,sindex(c,b,J[i].size(),J[j].size()),J[Jidx].size())] -= 2 * gamma * inv_qi;
            for (int d=b;d<J[j].size();d++) hess[hidx][hindex(sab,sindex(a,d,J[i].size(),J[j].size()),J[Jidx].size())] -= 2 * gamma * inv_qj;
        
        } }
        
        for (int k=0;k<hess[hidx].size();k++) hess[hidx][k] += 2 * gamma * inv_qi * inv_qj;
        
        // Add contribution of "zero" couplings to regularization term
        
        for (int a=0;a<J[i].size();a++) gamma_J_squared += pow((sumall * inv_qi * inv_qj) - (suma[a] * inv_qj), 2.0);
        for (int b=0;b<J[j].size();b++) gamma_J_squared += pow((sumall * inv_qi * inv_qj) - (sumb[b] * inv_qi), 2.0);
        gamma_J_squared += pow(sumall * inv_qi * inv_qj, 2.0);
        
    } }
	
	// Return the result
	
    gamma_J_squared *= gamma;
	S               += gamma_J_squared;
	
}


// Modifies the supplied model entropy given a set of couplings, with L2 regularization

void regularizeS_L2_gradOnly_GI(const Vector &J, double &S, Vector &grad, const Vector &p, double gamma, double gammah) {
    
    int length=(int) sizetolength(J.size());
    double gamma_J_squared=0;
    
    // Regularization for fields
    
    for (int i=0;i<length;i++) { for (int j=0;j<J[i].size();j++) {
            
        gamma_J_squared += gammah * J[i][j] * J[i][j];
        grad[i][j]      += 2 * gammah * gamma * J[i][j];
    
    } }

    // GAUGE INVARIANT Regularization for couplings
    
    for (int i=0;i<length;i++) { for (int j=i+1;j<length;j++) {
    
        // Get quantities for computing K_ij^ab
    
        int Jidx = index(i,j,length);
        double inv_qi = 1/((double)J[i].size()+1);
        double inv_qj = 1/((double)J[j].size()+1);
    
        double suma[J[i].size()];
        double sumb[J[j].size()];
        double sumall=0;
        
        for (int a=0;a<J[i].size();a++) suma[a]=0;
        for (int b=0;b<J[j].size();b++) sumb[b]=0;
        
        for (int a=0;a<J[i].size();a++) { for (int b=0;b<J[j].size();b++) {
            
            suma[a] += J[Jidx][sindex(a,b,J[i].size(),J[j].size())];
            sumb[b] += J[Jidx][sindex(a,b,J[i].size(),J[j].size())];
            
        } }
        
        for (int a=0;a<J[i].size();a++) sumall += suma[a];
        
        // Compute K_ij^ab and contribution to regularization term and Hessian
        
        for (int a=0;a<J[i].size();a++) { for (int b=0;b<J[j].size();b++) {
        
            int sab  = sindex(a,b,J[i].size(),J[j].size());
            double K = J[Jidx][sab] - (suma[a] * inv_qj) - (sumb[b] * inv_qi) + (sumall * inv_qi * inv_qj);
    
            gamma_J_squared += K * K;
            grad[Jidx][sab] += 2 * gamma * K;
        
        } }
        
        // Add contribution of "zero" couplings to regularization term
        
        for (int a=0;a<J[i].size();a++) gamma_J_squared += pow((sumall * inv_qi * inv_qj) - (suma[a] * inv_qj), 2.0);
        for (int b=0;b<J[j].size();b++) gamma_J_squared += pow((sumall * inv_qi * inv_qj) - (sumb[b] * inv_qi), 2.0);
        gamma_J_squared += pow(sumall * inv_qi * inv_qj, 2.0);
        
    } }

	// Return the result
    
    gamma_J_squared *= gamma;
	S               += gamma_J_squared;
	
}

