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
    
    // Test array-array eraseInPlace
    
    int v0 [1]      = {0};
    int len0        = 1;
    int result0 [1] = {0};
    
    eraseInPlace(v0,len0,0);
    assert(v0[0]==result0[0] && "tools - eraseInPlace, array gave an unexpected result");
    
    int v1 [4]      = {0, 1, 3, 4};
    int len1        = 4;
    int result1 [4] = {0, 3, 4, 4};
    
    eraseInPlace(v1,len1,1);
    assert(v1[0]==result1[0] && "tools - eraseInPlace, array gave an unexpected result");
    assert(v1[1]==result1[1] && "tools - eraseInPlace, array gave an unexpected result");
    assert(v1[2]==result1[2] && "tools - eraseInPlace, array gave an unexpected result");
    assert(v1[3]==result1[3] && "tools - eraseInPlace, array gave an unexpected result");
    v1[1] = 1;
    v1[2] = 3;
    v1[3] = 4;
    
    // Test vector eraseInPlace
    
    std::vector<int> vv0, vresult0;
    for (int i=0;i<len0;i++)   { vv0.push_back(v0[i]);               }
    
    eraseInPlace(vv0,0);
    assert(vv0==vresult0 && "tools - eraseInPlace, vector gave an unexpected result");
    
    int result1temp [3] = {0, 3, 4};
    std::vector<int> vv1(len1,0), vresult1(len1-1,0);
    for (int i=0;i<len1;i++)   { vv1[i]=v1[i];               }
    for (int i=0;i<len1-1;i++) { vresult1[i]=result1temp[i]; }
    
    eraseInPlace(vv1,1);
    assert(vv1==vresult1 && "tools - eraseInPlace, vector gave an unexpected result");
    
	return 0;

}
