//
//  matrix_unittest.cpp
//  MLib
//
//  Created by Iaroslav Omelianenko on 1/22/15.
//  Copyright (c) 2015 yaric. All rights reserved.
//
#include "unit++/unit++.h"
#include <sstream>
#include <vector>

#include "defines.h"
#include "matrix.h"
#include "utils.h"
#include "dim_reduction_methods.h"
#include "features_normalization.h"
#include "io.h"
#include "gradient_boosting_tree.h"
#include "RF_regression.h"
#include "estimators.h"
#include "ridge_regression.h"
#include "random.h"

using namespace nologin;
using namespace nologin::math;
using namespace nologin::utils;
using namespace nologin::preprocessing;
using namespace nologin::tree;

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
    
    void checkExpandMatrix() {
        Printf("Check Matrix Expand Lineary+++++++++++++++++++++++++++\n");
        Matrix A(2, {1, 2,
                    3, 4,
                    5, 6,
                    7, 8});
        print(A);
        
        Matrix B = expandMatrixLineary(A);
        print(B);
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
        
        add("checkExpandMatrix", unitpp::testcase(this, "Check Matrix Expand Lineary", &TestMatrix::checkExpandMatrix));
        
        // add this suite to the main suite
        suite::main().add("Matrix", this);
    }
};

class RandomTest : public unitpp::suite {
    
    void checkNextDouble() {
        Printf("Check next double random number+++++++++++++++++++++++++++\n");
        RNG rnd;
        for (int i = 0; i < 1000; i++) {
            double d = rnd.nextDouble();
            Assert(d < 1 && d >= 0, "Random value if out of scope: %f", d);
            Printf("%f\n", d);
        }
    }
    
    
public:
    RandomTest() : suite("Random test") {
        add("checkNextDouble", unitpp::testcase(this, "Check next double random number", &RandomTest::checkNextDouble));
        
        // add this suite to the main suite
        suite::main().add("Random", this);
    }
};

class TestGBTModelSerialization : public unitpp::suite {
    void checkNodeSerializationTest() {
        Printf("The Node tree serialization test+++++++++++++++++++++++++++\n");
        Node *n1 = createNode();
        
        ostringstream so;
        saveNode(n1, so);
        
        string output = so.str();
        Printf("Serialized:\t\t%s\n", output.c_str());
        
        // deserialize
        istringstream si(output);
        Node *n_out = loadNode(si);
        unitpp::assert_eq("Root Node value", n1->m_node_value, n_out->m_node_value);
        unitpp::assert_true("Root Node left value", abs(n1->m_left_child->m_node_value - n_out->m_left_child->m_node_value) < 0.0001);
        unitpp::assert_true("Root Node right value", abs(n1->m_right_child->m_node_value - n_out->m_right_child->m_node_value) < 0.0001);
        
        ostringstream os;
        saveNode(n_out, os);
        Printf("Deserialized:\t%s\n", os.str().c_str());
        
        delete n1;
    }
    
    void checkSaveRegressionTree() {
        Printf("The regression tree serialization test+++++++++++++++++++++++++++\n");
        Node *root = createNode();
        RegressionTree tree;
        tree.m_root = root;
        tree.m_current_depth = 11;
        tree.m_max_depth = 20;
        tree.m_min_nodes = 5;
        
        
        // serialize
        ostringstream so;
        saveRegressionTree(tree, so);
        string output = so.str();
        Printf("Serialized:\t\t%s\n", output.c_str());
        
        // deserialize
        istringstream si(output);
        RegressionTree d_tree;
        loadRegressionTree(d_tree, si);
        unitpp::assert_eq("Tree depth", d_tree.m_current_depth, tree.m_current_depth);
        unitpp::assert_eq("Tree max depth", d_tree.m_max_depth, tree.m_max_depth);
        unitpp::assert_eq("Tree min nodes", d_tree.m_min_nodes, tree.m_min_nodes);
        unitpp::assert_eq("Tree Root Node value", d_tree.m_root->m_node_value, tree.m_root->m_node_value);
        
        ostringstream os;
        saveRegressionTree(d_tree, os);
        Printf("Deserialized:\t%s\n", os.str().c_str());
    }
    
    void checkGPFSerializationTest() {
        Printf("The Gradient Prediction Forest serialization test+++++++++++++++++++++++++++\n");
        Node *root = createNode();
        RegressionTree tree;
        tree.m_root = root;
        tree.m_current_depth = 11;
        tree.m_max_depth = 20;
        tree.m_min_nodes = 5;
        
        PredictionForest forest(0.5);
        forest.m_init_value = 10;
        forest.m_trees = {tree};
        
        // serialize
        ostringstream so;
        savePredictionForest(forest, so);
        string output = so.str();
        Printf("Serialized:\t\t%s\n", output.c_str());
        
        // deserialize
        istringstream si(output);
        PredictionForest d_forest(0.01);
        loadPredictionForest(d_forest, si);
        unitpp::assert_eq("Init value", d_forest.m_init_value, forest.m_init_value);
        unitpp::assert_eq("Learning rate", d_forest.m_combine_weight, forest.m_combine_weight);
        unitpp::assert_eq("Trees count", d_forest.m_trees.size(), forest.m_trees.size());
        
        RegressionTree d_tree = d_forest.m_trees[0];
        unitpp::assert_eq("Tree depth", d_tree.m_current_depth, tree.m_current_depth);
        unitpp::assert_eq("Tree max depth", d_tree.m_max_depth, tree.m_max_depth);
        unitpp::assert_eq("Tree min nodes", d_tree.m_min_nodes, tree.m_min_nodes);
        unitpp::assert_eq("Tree Root Node value", d_tree.m_root->m_node_value, tree.m_root->m_node_value);
        
        ostringstream os;
        savePredictionForest(d_forest, os);
        Printf("Deserialized:\t%s\n", os.str().c_str());
    }
    
    
    Node *createNode() {
        Node *n5 = new Node(.5, 5, .5, .5);
        Node *n4 = new Node(.4, 4, .4, .4);
        Node *n3 = new Node(.3, 3, .3, .3);
        n3->m_left_child = n4;
        n3->m_right_child = n5;
        Node *n2 = new Node(.2, 2, .2, .2);
        Node *n1 = new Node(.1, 1, .1, .1);
        n1->m_left_child = n3;
        n1->m_right_child = n2;
        return n1;
    }
    
public:
    TestGBTModelSerialization() : suite("The Gradient Boosting Tree Model Serialization Test Suite") {
        add("checkNodeSerializationTest", unitpp::testcase(this, "The Node tree serialization test", &TestGBTModelSerialization::checkNodeSerializationTest));
        add("checkSaveRegressionTree", unitpp::testcase(this, "The regression tree serialization test", &TestGBTModelSerialization::checkSaveRegressionTree));
        add("checkGPFSerializationTest", unitpp::testcase(this, "The Gradient Prediction Forest serialization test", &TestGBTModelSerialization::checkGPFSerializationTest));
        
        // add this suite to the main suite
        suite::main().add("The Model Serialization", this);
    }
};

RandomTest *test = new RandomTest();
TestMatrix* theTest = new TestMatrix();
TestGBTModelSerialization *modelTest = new TestGBTModelSerialization();
