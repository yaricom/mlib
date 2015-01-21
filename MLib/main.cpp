//
//  main.cpp
//  MLib
//
//  Created by Iaroslav Omelianenko on 1/21/15.
//  Copyright (c) 2015 yaric. All rights reserved.
//

#include <iostream>

#include "defines.h"
#include "matrix.h"
#include "utils.h"

using namespace std;

void checkMatrixMean() {
    VD data = {0, 1, 1, 2, 3, 2, 1, 3, 2, 4, 2, 2};
    Matrix m(3, data);
    print(m);
    Matrix mean = m.mean();
    print(mean);
    
    VD test = {1.75,	2.25,	1.75};
    Matrix tm(3, test);
    Assert(mean == tm, "Calculated matrix mean is wrong!");
}

int main(int argc, const char * argv[]) {
    // check matrix mean
    checkMatrixMean();
    
    
    
    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}
