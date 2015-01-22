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
#include "dim_reduction_methods.h"

using namespace std;

void checkMatrixMean() {
    Printf("Matrix mean+++++++++++++++++++++++++++\n");
    VD data = {0, 1, 1, 2, 3, 2, 1, 3, 2, 4, 2, 2};
    Matrix m(3, data);
    print(m);
    Matrix mean = m.mean();
    print(mean);
    
    VD test = {1.75,	2.25,	1.75};
    Matrix tm(3, test);
    Assert(mean == tm, "Calculated matrix mean is wrong!");
}

void checkMatrixStd() {
    Printf("Matrix STD+++++++++++++++++++++++++++\n");
    VD data = {1, 5, 9, 7, 15, 22};
    Matrix m(3, data);
    print(m);
    
    Matrix std = m.std();
    print(std);
    
    VD test = {4.2426406871192848, 7.0710678118654755, 9.1923881554251174};
    Matrix tm(3, test);
    Assert(std == tm, "Calculated matrix std deviation is wrong!");
}

void checkRawAccess() {
    Printf("Raw access+++++++++++++++++++++++++++\n");
    VD data = {0, 1, 1, 2, 3, 2, 1, 3, 2, 4, 2, 2};
    Matrix m(3, data);
    print(m);
    
    double expected, test;
    size_t rows = data.size() / 3;
    // check read
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < 3; j++) {
            expected = data[i * 3 + j];
            test = m(i, j);
            Assert(test == expected, "Expected: %f but was: %f", expected, test);
        }
    }
    
    // check write
    int val = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < 3; j++) {
            m(i, j) = val++;
        }
    }
    VD testVector;
    m.rowPackedCopy(testVector);
    for (int i = 0; i < testVector.size(); i++) {
        test = testVector[i];
        Assert(test == i, "Expected: %f but was: %f", i, test);
    }
}

void checkCopy() {
    Printf("Check copy+++++++++++++++++++++++++++\n");
    VD data = {0, 1, 1, 2, 3, 2, 1, 3, 2, 4, 2, 2};
    Matrix m(3, data);
    print(m);
}

int main(int argc, const char * argv[]) {
    // check matrix mean
    checkMatrixMean();
    // check standard deviation
    checkMatrixStd();
    // check raw access to the matrix data
    checkRawAccess();
    
    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}
