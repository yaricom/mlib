//
//  matrix_unittest.cpp
//  MLib
//
//  Created by Iaroslav Omelianenko on 1/22/15.
//  Copyright (c) 2015 yaric. All rights reserved.
//

#include "unit++/unit++.h"
#include "defines.h"
#include "matrix.h"
#include "utils.h"

class TestMatrix : public unitpp::suite {
    
    void checkMatrixMean() {
        Printf("Matrix mean+++++++++++++++++++++++++++\n");
        VD data = {0, 1, 1, 2, 3, 2, 1, 3, 2, 4, 2, 2};
        Matrix m(3, data);
        print(m);
        Matrix mean = m.mean();
        print(mean);
        
        VD test = {1.75,	2.25,	1.75};
        Matrix tm(3, test);
        unitpp::assert_true("Calculated matrix mean is wrong!", mean == tm);
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
        unitpp::assert_true("Calculated matrix std deviation is wrong!", std == tm);
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
                unitpp::assert_eq(spf("Expected value not equal to readen at row: %i, col: %i", i, j), test, expected);
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
            unitpp::assert_eq(spf("Value not expected"), test, i);
        }
    }
    
    void checkCopy() {
        Printf("Check copy+++++++++++++++++++++++++++\n");
        VD data = {0, 1, 1, 2, 3, 2, 1, 3, 2, 4, 2, 2};
        Matrix m(3, data);
        print(m);
#warning implement this!
    }
    
    void checkScaleMinMax() {
        VD data = {1, 5, 9, 7, 15, 22};
        Matrix m(3, data);
        print(m);
#warning implement this!
    }
    
public:
    TestMatrix() : suite("The Matrix Test Suite") {
        
        add("checkMatrixMean", unitpp::testcase(this, "Matrix Mean test", &TestMatrix::checkMatrixMean));
        add("checkMatrixStd", unitpp::testcase(this, "Matrix STD test", &TestMatrix::checkMatrixStd));
        add("checkRawAccess", unitpp::testcase(this, "Matrix Raw access test", &TestMatrix::checkRawAccess));
        add("checkCopy", unitpp::testcase(this, "Matrix Check copy test", &TestMatrix::checkCopy));
        
        // add this suite to the main suite
        suite::main().add("Matrix", this);
    }
};
TestMatrix* theTest = new TestMatrix();