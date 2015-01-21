//
//  dim_reduction_methods.h
//  activemoleculesC++
//
//  Created by Iaroslav Omelianenko on 1/21/15.
//  Copyright (c) 2015 yaric. All rights reserved.
//

#ifndef activemoleculesC___dim_reduction_methods_h
#define activemoleculesC___dim_reduction_methods_h

#include "Definitions.h"
#include "Matrix.h"
#include "Random.h"
#include "random_projection.h"

inline Matrix embedRandomProjection(const Matrix &data, const IndexType p_target_dimension) {
    IndexType current_dimension = data.getColumnDimension();
    
    Matrix projection_matrix = gaussian_projection_matrix(current_dimension, p_target_dimension);
    DenseVector mean_vector = compute_mean(begin, end, features, current_dimension);
    
    tapkee::ProjectingFunction projecting_function(new tapkee::MatrixProjectionImplementation(projection_matrix,mean_vector));
    return TapkeeOutput(project(projection_matrix,mean_vector,begin,end,features,current_dimension), projecting_function);
}

#endif
