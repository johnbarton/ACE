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
    
    // Test array-array L1 norm
    
    double v1 [5] = {0, -1, 2, -3, 4};
    double result = 10;
    int len       = 5;
    
    assert(L1(v1,len)==result && "tools - L1, array gave an unexpected result");
    
    // Test vector L1 norm
    
    std::vector<double> vv1;
    for (int i=0;i<len;i++) { vv1.push_back(v1[i]); }
    
    assert(L1(vv1)==result && "tools - L1, vector gave an unexpected result");
    
    // Test Vector L1 norm
    
    Vector vvv1;
    for (int i=0;i<len;i++) { vvv1.push_back(vv1); }
    
    assert(L1(vvv1)==(len*result) && "tools - L1, Vector gave an unexpected result");
    
	return 0;

}
