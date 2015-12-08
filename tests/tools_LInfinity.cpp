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
    
    // Test array L Infinity norm
    
    double v1 [5] = {0, -1, 2, 3, -4};
    double v2 [5] = {4, -3, -2, 1, 0};
    double result = 4;
    int len       = 5;
    
    assert(LInfinity(v1,len)==result && "tools - LInfinity, array gave an unexpected result");
    assert(LInfinity(v2,len)==result && "tools - LInfinity, array gave an unexpected result");
    
    // Test vector L Infinity norm
    
    std::vector<double> vv1, vv2;
    for (int i=0;i<len;i++) { vv1.push_back(v1[i]); vv2.push_back(v2[i]); }
    
    assert(LInfinity(vv1)==result && "tools - LInfinity, vector gave an unexpected result");
    assert(LInfinity(vv2)==result && "tools - LInfinity, vector gave an unexpected result");
    
    // Test Vector L Infinity norm
    
    Vector vvv1, vvv2;
    for (int i=0;i<len;i++) { vvv1.push_back(vv1); vvv2.push_back(vv2); }
    
    assert(LInfinity(vvv1)==result && "tools - LInfinity, Vector gave an unexpected result");
    assert(LInfinity(vvv2)==result && "tools - LInfinity, Vector gave an unexpected result");
    
    // Test sparse Vector L Infinity norm
    
    IntVector sparse;
    for (int i=0;i<len;i++) { sparse.push_back(std::vector<int>()); sparse[i].push_back(i); }
    
    assert(LInfinity(vvv1,sparse)==result && "tools - LInfinity, sparse Vector gave an unexpected result");
    assert(LInfinity(vvv2,sparse)==result && "tools - LInfinity, sparse Vector gave an unexpected result");
    
	return 0;

}
