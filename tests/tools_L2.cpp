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

#include "../src/tools.h"          // Numerical tools

int main(int argc, char *argv[]) {
    
    // Test array-array L2 norm
    
    double v1 [5]      = {0, 3, 3, 3, 3};
    double v1short [3] = {0, 3, 3};
    double result      = 6;
    int len            = 5;
    int lenshort       = 3;
    
    assert(L2(v1,len)==result && "tools - L2, array gave an unexpected result");
    
    // Test vector L2 norm
    
    std::vector<double> vv1;
    for (int i=0;i<len;i++) { vv1.push_back(v1[i]); }
    
    assert(L2(vv1)==result && "tools - L2, vector gave an unexpected result");
    
    // Test Vector L2 norm
    
    std::vector<double> vv1short;
    for (int i=0;i<lenshort;i++) { vv1short.push_back(v1short[i]); }
    
    Vector vvv1;
    for (int i=0;i<2;i++) { vvv1.push_back(vv1short); }
    
    assert(L2(vvv1)==result && "tools - L2, Vector gave an unexpected result");
    
	return 0;

}
