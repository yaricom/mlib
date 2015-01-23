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
        Vector vmean = mean(m);
        print(vmean);
        
        VD test = {1.75, 2.25, 1.75};
        Vector tm(test);
        unitpp::assert_true("Calculated matrix mean is wrong!", vmean == tm);
    }
    
    void checkMatrixVariance() {
        Printf("Matrix Variance+++++++++++++++++++++++++++\n");
        VD data = {4, -2, 1, 9, 5, 7};
        Matrix m(3, data);
        print(m);
        
        Vector var = variance(m);
        print(var);
        
        VD test = {12.5000, 24.5000, 18.0000};
        Vector tm(test);
        unitpp::assert_true("Calculated matrix variance is wrong!", var == tm);
    }
    
    void checkMatrixStd() {
        Printf("Matrix STD+++++++++++++++++++++++++++\n");
        VD data = {1, 5, 9, 7, 15, 22};
        Matrix m(3, data);
        print(m);
        
        Vector std = stdev(m);
        print(std);
        
        VD test = {4.2426406871192848, 7.0710678118654755, 9.1923881554251174};
        Vector tm(test);
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
            unitpp::assert_eq(spf("Value found: %i, but was expected: %i", test, i), test, i);
        }
    }
    
    void checkCopy() {
        Printf("Check copy+++++++++++++++++++++++++++\n");
        VD data = {0, 1, 1, 2, 3, 2, 1, 3, 2, 4, 2, 2};
        Matrix m(3, data);
        print(m);
        
        Matrix test = m;
        print(test);
        
        unitpp::assert_true("Copied matrix has wrong content!", test == m);
    }
    
    void checkScaleMinMax() {
        Printf("Check scale Min-Max+++++++++++++++++++++++++++\n");
        VD data = {1, 5, 9, 2, 3, 39, 7, 15, 22};
        Matrix m(3, data);
        print(m);
        
        double min = -1;
        double max = 1;
        m.scaleMinMax(min, max);
        print(m);
        
        double tMin = 1000, tMax = -1000;
        for (int i = 0 ; i < m.rows(); i++) {
            for (int j = 0; j < m.cols(); j++) {
                if (tMax < m(i, j)) {
                    tMax = m(i, j);
                }
                if (tMin > m(i, j)) {
                    tMin = m(i, j);
                }
            }
        }
        unitpp::assert_true(spf("Expected matrix min should be greater than: %f, but was: %f", min, tMin), tMin >= min);
        unitpp::assert_true(spf("Expected matrix max should be less than: %f, but was: %f", max, tMax), tMax <= max);
    }
    
    void checkMatrixProduct() {
        Printf("Check Matrix Product+++++++++++++++++++++++++++\n");
        VD dataA = {1, 2, 3, 4, 5, 6, 7, 8};
        Matrix A(2, dataA);
        print(A);
        
        VD dataB = {1, 2, 3, 4, 5, 6};
        Matrix B(3, dataB);
        print(B);
        
        Matrix C = A.matmul(B);
        print(C);
        
        VD dataCheck = {9, 12, 15, 19, 26, 33, 29, 40, 51, 39, 54, 69};
        Matrix test(3, dataCheck);
        
        unitpp::assert_true("Calculated matrix product is wrong!", C == test);
    }
    
public:
    TestMatrix() : suite("The Matrix Test Suite") {
        
        add("checkMatrixMean", unitpp::testcase(this, "Matrix Mean test", &TestMatrix::checkMatrixMean));
        add("checkMatrixVariance", unitpp::testcase(this, "Matrix Variance test", &TestMatrix::checkMatrixVariance));
        add("checkMatrixStd", unitpp::testcase(this, "Matrix STD test", &TestMatrix::checkMatrixStd));
        add("checkRawAccess", unitpp::testcase(this, "Matrix Raw access test", &TestMatrix::checkRawAccess));
        add("checkCopy", unitpp::testcase(this, "Matrix Check copy test", &TestMatrix::checkCopy));
        add("checkScaleMinMax", unitpp::testcase(this, "Matrix Scale Min-Max test", &TestMatrix::checkScaleMinMax));
        add("checkMatrixProduct", unitpp::testcase(this, "Matrix Product test", &TestMatrix::checkMatrixProduct));
        
        // add this suite to the main suite
        suite::main().add("Matrix", this);
    }
};
TestMatrix* theTest = new TestMatrix();