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

/**
 * Makes Gaussian Random Projection of provided matrix to the target dimension, reducing features dimension.
 */
Matrix& randomProjection(const Matrix &data, const int target_dimension) {
    size_t current_dimension = data.cols();
    
    Matrix projection_matrix = Matrix::gaussianRandom((int)current_dimension, target_dimension);
    
    Matrix &projected = data.matmul(projection_matrix);
    
    return projected;
}

#endif
