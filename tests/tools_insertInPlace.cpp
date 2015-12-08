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
    
    // Test array-array insertInPlace
    
    int v0 [1]      = {0};
    int len0        = 1;
    int result0 [1] = {1};
    
    insertInPlace(v0,len0,1);
    assert(v0[0]==result0[0] && "tools - insertInPlace, array gave an unexpected result");
    v0[0] = 0;
    
    int v1 [4]      = {0, 1, 3, 4};
    int len1        = 4;
    int result1 [4] = {0, 1, 2, 3};
    
    insertInPlace(v1,len1,2);
    assert(v1[0]==result1[0] && "tools - insertInPlace, array gave an unexpected result");
    assert(v1[1]==result1[1] && "tools - insertInPlace, array gave an unexpected result");
    assert(v1[2]==result1[2] && "tools - insertInPlace, array gave an unexpected result");
    assert(v1[3]==result1[3] && "tools - insertInPlace, array gave an unexpected result");
    v1[2] = 3;
    v1[3] = 4;
    
    // Test vector insertInPlace
    
    int result0temp [2] = {0, 1};
    std::vector<int> vv0, vresult0;
    for (int i=0;i<len0;i++)   { vv0.push_back(v0[i]);               }
    for (int i=0;i<len0+1;i++) { vresult0.push_back(result0temp[i]); }
    
    insertInPlace(vv0,1);
    assert(vv0==vresult0 && "tools - insertInPlace, vector gave an unexpected result");
    
    int result1temp [5] = {0, 1, 2, 3, 4};
    std::vector<int> vv1(len1,0), vresult1(len1+1,0);
    for (int i=0;i<len1;i++)   { vv1[i]=v1[i];               }
    for (int i=0;i<len1+1;i++) { vresult1[i]=result1temp[i]; }

    insertInPlace(vv1,2);
    //printf("test: ");
    //for (int i=0;i<vv1.size();i++) printf("%d ",vv1[i]);
    printf("test: ");
    for (int i=0;i<len1;i++) printf("%d ",v1[i]);
    printf("\texpected: ");
    for (int i=0;i<vresult1.size();i++) printf("%d ",vresult1[i]);
    assert(vv1==vresult1 && "tools - insertInPlace, vector gave an unexpected result");
    
	return 0;

}
