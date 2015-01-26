//
//  dim_reduction_methods.h
//  activemoleculesC++
//
//  Created by Iaroslav Omelianenko on 1/21/15.
//  Copyright (c) 2015 yaric. All rights reserved.
//

#ifndef activemoleculesC___dim_reduction_methods_h
#define activemoleculesC___dim_reduction_methods_h

#include "matrix.h"
#include "random.h"


Matrix gaussian_projection_matrix(int target_dimension, int current_dimension) {
    Matrix projection_matrix(target_dimension, current_dimension);
    
    double div = sqrt(static_cast<double>(target_dimension));
    for (int i = 0; i < target_dimension; ++i) {
        for (int j = 0; j < current_dimension; ++j) {
            projection_matrix(i, j) = gaussian_random() / div;
        }
    }
    
    return projection_matrix;
}

/**
 * Makes Gaussian Random Projection of provided matrix to the target dimension, reducing features dimension.
 */
Matrix& randomProjection(const Matrix &data, const int target_dimension) {
    size_t current_dimension = data.cols();
    
    Matrix projection_matrix = gaussian_projection_matrix((int)current_dimension, target_dimension);
    
    Matrix &projected = data.matmul(projection_matrix);
    
    return projected;
}

#endif
