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
#include "../utils/random.h"

inline Matrix gaussian_projection_matrix(IndexType target_dimension, IndexType current_dimension) {
    Matrix projection_matrix(target_dimension,current_dimension);
    
    for (IndexType i=0; i<target_dimension; ++i) {
        for (IndexType j=0; j<current_dimension; ++j) {
            projection_matrix.set(i, j, gaussian_random()/sqrt(static_cast<ScalarType>(target_dimension)));
        }
    }
    
    return projection_matrix;
}

#endif
