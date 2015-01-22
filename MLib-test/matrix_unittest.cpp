//
//  matrix_unittest.cpp
//  MLib
//
//  Created by Iaroslav Omelianenko on 1/22/15.
//  Copyright (c) 2015 yaric. All rights reserved.
//

#include <stdio.h>
#include <gtest/gtest.h>

#include "matrix.h"
#include "utils.h"

// TEST has two parameters: the test case name and the test name.
// After using the macro, you should define your test logic between a
// pair of braces.  You can use a bunch of macros to indicate the
// success or failure of a test.  EXPECT_TRUE and EXPECT_EQ are
// examples of such macros.  For a complete list, see gtest.h.

// Tests matrix mean calculation.
TEST(Matrix, MatrixMean) {
    Printf("Matrix mean+++++++++++++++++++++++++++\n");
    VD data = {0, 1, 1, 2, 3, 2, 1, 3, 2, 4, 2, 2};
    Matrix m(3, data);
    print(m);
    Matrix mean = m.mean();
    print(mean);
    
    VD test = {1.75,	2.25,	1.75};
    Matrix tm(3, test);
    EXPECT_EQ(mean, tm) << "Calculated matrix mean is wrong!";
}