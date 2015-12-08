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
    
    // Test array-array inner product
    
    double v1 [5] = {0, 1, 2, 3, 4};
    double v2 [5] = {4, 3, 2, 1, 0};
    double result = 10;
    int len       = 5;
    
    assert(innerProduct(v1,v2,len)==result && "tools - innerProduct, array-array gave an unexpected result");
    
    // Test vector-vector inner product
    
    std::vector<double> vv1, vv2;
    for (int i=0;i<len;i++) { vv1.push_back(v1[i]); vv2.push_back(v2[i]); }
    
    assert(innerProduct(vv1,vv2)==result && "tools - innerProduct, vector-vector gave an unexpected result");
    
    // Test Vector-Vector inner product
    
    Vector vvv1, vvv2;
    for (int i=0;i<len;i++) { vvv1.push_back(vv1); vvv2.push_back(vv2); }
    
    assert(innerProduct(vvv1,vvv2)==(len*result) && "tools - innerProduct, Vector-Vector gave an unexpected result");
    
    // Test sparse Vector-Vector inner product
    
    IntVector sparse;
    for (int i=0;i<len;i++) { sparse.push_back(std::vector<int>()); sparse[i].push_back(i); }
    
    assert(innerProduct(vvv1,vvv2,sparse)==result && "tools - innerProduct, sparse Vector-Vector gave an unexpected result");
    
	return 0;

}
