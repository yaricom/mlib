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
#include "dim_reduction_methods.h"
#include "features_normalization.h"
#include "io.h"

using namespace nologin;
using namespace nologin::math;
using namespace nologin::utils;
using namespace nologin::preprocessing;

class TestMatrix : public unitpp::suite {
    
    void checkMatrixMean() {
        Printf("Matrix mean+++++++++++++++++++++++++++\n");
        Matrix m(3, {0, 1, 1,
                    2, 3, 2,
                    1, 3, 2,
                    4, 2, 2});
        print(m);
        Vector vmean = m.mean();
        print(vmean);
        
        Vector tm({1.75, 2.25, 1.75});
        unitpp::assert_true("Calculated matrix mean is wrong!", vmean == tm);
    }
    
    void checkMatrixVariance() {
        Printf("Matrix Variance+++++++++++++++++++++++++++\n");
        Matrix m(3, {4, -2, 1,
                    9, 5, 7});
        print(m);
        
        Vector var = m.variance(1);
        print(var);
        
        Vector tm({12.5000, 24.5000, 18.0000});
        unitpp::assert_true("Calculated matrix variance is wrong!", var == tm);
    }
    
    void checkMatrixStd() {
        Printf("Matrix STD+++++++++++++++++++++++++++\n");
        Matrix m(3, {1, 5, 9,
                    7, 15, 22});
        print(m);
        
        Vector stdM = m.stdev(1);
        print(stdM);
        
        Vector tm({4.242640, 7.071067, 9.192388});
        unitpp::assert_true("Calculated matrix std deviation is wrong!", stdM.similar(tm, 0.000001));
    }
    
    void checkRawAccess() {
        Printf("Raw access+++++++++++++++++++++++++++\n");
        VD data = { 0, 1, 1,
                    2, 3, 2,
                    1, 3, 2,
                    4, 2, 2};
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
        Matrix m(3, {0, 1, 1,
                    2, 3, 2,
                    1, 3, 2,
                    4, 2, 2});
        print(m);
        
        Matrix test = m;
        print(test);
        
        unitpp::assert_true("Copied matrix has wrong content!", test == m);
    }
    
    void checkScaleMinMax() {
        Printf("Check scale Min-Max+++++++++++++++++++++++++++\n");
        Matrix m(3, {1, 5, 9,
                    2, 3, 39,
                    7, 15, 22});
        print(m);
        
        double min = -1;
        double max = 1;
        Matrix testM = m;
        scaleMinMax(min, max, testM);
        print(testM);
        
        double tMin = 1000, tMax = -1000;
        for (int i = 0 ; i < testM.rows(); i++) {
            for (int j = 0; j < testM.cols(); j++) {
                if (tMax < testM(i, j)) {
                    tMax = testM(i, j);
                }
                if (tMin > testM(i, j)) {
                    tMin = testM(i, j);
                }
            }
        }
        unitpp::assert_true(spf("Expected matrix min should be greater than: %f, but was: %f", min, tMin), tMin >= min);
        unitpp::assert_true(spf("Expected matrix max should be less than: %f, but was: %f", max, tMax), tMax <= max);
        
        // check scale by indices
        VI indices = {0, 2};
        Matrix testM2 = m;
        scaleMinMax(min, max, indices, testM2);
        print(testM2);
        Matrix checkM(3, {-1, 5, -1,
                        -0.666666, 3, 1,
                        1, 15, -0.133333});
        
        unitpp::assert_true("Scale by min/max with column indices failed", testM2.similar(checkM, 0.000001));
    }
    
    void checkStdScale() {
        Printf("Check standard scale+++++++++++++++++++++++++++\n");
        Matrix m(3, {1, -2, 2,
                    2, 0, 0,
                    0, 1, -1});
        print(m);
        
        Matrix testM = m;
        stdScale(testM);
        print(testM);
        
        Vector stdev = testM.stdev(1);
        Vector checkDev({1, 1, 1});

        unitpp::assert_true("Wrong std deviance values after matrix scale", checkDev.similar(stdev, .000001));
        
        Vector stdMean = testM.mean();
        Vector checkMean({0, 0, 0});
        unitpp::assert_true("Wrong mean values after matrix scale", stdMean.similar(checkMean, .000001));
        
        // check scale by indices
        VI indices = {0, 2};
        Matrix testM2 = m;
        stdScale(indices, testM2);
        print(testM2);
        
        Vector stdev2 = testM2.stdev(1);
        unitpp::assert_true("Wrong std deviance values after matrix scale by indices", abs(stdev2[0] - stdev2[2]) <  .000001);
        Vector stdMean2 = testM2.mean();
        unitpp::assert_true("Wrong mean values after matrix scale by indices", stdMean2[0] == stdMean2[2] && stdMean2[0] == 0);
    }
    
    void checkMatrixProduct() {
        Printf("Check Matrix Product+++++++++++++++++++++++++++\n");
        Matrix A(2, {1, 2,
                    3, 4,
                    5, 6,
                    7, 8});
        print(A);
        
        Matrix B(3, {1, 2, 3,
                    4, 5, 6});
        print(B);
        
        Matrix C = A.matmul(B);
        print(C);
        
        Matrix test(3, {9, 12, 15,
                        19, 26, 33,
                        29, 40, 51,
                        39, 54, 69});
        
        unitpp::assert_true("Calculated matrix product is wrong!", C == test);
    }
    
    void checkGaussianRandomProjection() {
        Printf("Check Gaussian Random Projection+++++++++++++++++++++++++++\n");
        Matrix mat(4, {9, 12, 15, 21,
                        19, 26, 33, 35,
                        29, 40, 51, 14,
                        39, 54, 69, 67});
        print(mat);
        
        int target_dim = 3;
        Matrix res = randomProjection(mat, target_dim);
        print(res);
    }
    
    void checkCorrectOutliers() {
        Printf("Check Correct Outliers+++++++++++++++++++++++++++\n");
        Matrix mat(4, {-39, 12, 15, 21,
                        15, 16, 23, 25,
                        19, 26, 33, 35,
                         9,  7, 15, 22,
                        29, 40, 51, 14,
                        34, 20, 40, 30,
                        39, 54, 69, 167});
        print(mat);
        VI indices = {0, 3};
        correctOutliers(indices, mat);
        print(mat);
        
        unitpp::assert_true("Outliers was not corrected!", mat[0][0] > -99 && mat[3][3] < 167);
    }
    
    void checkMatrixStore() {
        Printf("Check Matrix Store/Load LibSVM+++++++++++++++++++++++++++\n");
        Matrix mat(4, {-39, 12, 0, 21,
                        15, 0, 23, 0,
                        0, 26, 33, 35,
                        9,  7, 0, 22,
                        29, 40, 51, 14,
                        34, 0, 0, 30,
                        39, 54, 69, 167});
        print(mat);
        
        string fileName = "/tmp/matrix.txt";
        storeMatrixAsLibSVM(fileName.c_str(), mat);
        
        // load matrix and check result
        Matrix &check = loadMatrixFromLibSVM(fileName.c_str());
        print(check);
        
        unitpp::assert_true("Failed to store/load matrix in LibSVM format", mat == check);
    }
    
public:
    TestMatrix() : suite("The Matrix Test Suite") {
        
        add("checkMatrixMean", unitpp::testcase(this, "Matrix Mean test", &TestMatrix::checkMatrixMean));
        add("checkMatrixVariance", unitpp::testcase(this, "Matrix Variance test", &TestMatrix::checkMatrixVariance));
        add("checkMatrixStd", unitpp::testcase(this, "Matrix STD test", &TestMatrix::checkMatrixStd));
        add("checkRawAccess", unitpp::testcase(this, "Matrix Raw access test", &TestMatrix::checkRawAccess));
        add("checkCopy", unitpp::testcase(this, "Matrix Check copy test", &TestMatrix::checkCopy));
        add("checkScaleMinMax", unitpp::testcase(this, "Matrix Scale Min-Max test", &TestMatrix::checkScaleMinMax));
        add("checkStdScale", unitpp::testcase(this, "Matrix Standard Scale test", &TestMatrix::checkStdScale));
        add("checkMatrixProduct", unitpp::testcase(this, "Matrix Product test", &TestMatrix::checkMatrixProduct));
        
        add("checkGaussianRandomProjection", unitpp::testcase(this, "Matrix Gaussian Random Projection test", &TestMatrix::checkGaussianRandomProjection));
        add("checkCorrectOutliers", unitpp::testcase(this, "Matrix Correct Outliers test", &TestMatrix::checkCorrectOutliers));
        
        add("checkMatrixStore", unitpp::testcase(this, "Matrix Store/Load LibSVM test", &TestMatrix::checkMatrixStore));
        
        // add this suite to the main suite
        suite::main().add("Matrix", this);
    }
};
TestMatrix* theTest = new TestMatrix();