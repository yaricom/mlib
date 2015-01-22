//
//  random_projection.h
//  activemoleculesC++
//
//  Created by Iaroslav Omelianenko on 1/21/15.
//  Copyright (c) 2015 yaric. All rights reserved.
//

#ifndef activemoleculesC___random_projection_h
#define activemoleculesC___random_projection_h

#include "Matrix.h"
#include "random.h"

inline Matrix gaussian_projection_matrix(int target_dimension, int current_dimension) {
    Matrix projection_matrix(target_dimension,current_dimension);
    
    double div = sqrt(static_cast<double>(target_dimension));
    for (int i = 0; i < target_dimension; ++i) {
        for (int j = 0; j < current_dimension; ++j) {
            projection_matrix(i, j) = gaussian_random() / div;
        }
    }
    
    return projection_matrix;
}

#endif
