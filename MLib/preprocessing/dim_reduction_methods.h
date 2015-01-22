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
#include "random_projection.h"

//template <class RandomAccessIterator, class FeatureVectorCallback>
//DenseMatrix project(const DenseMatrix& projection_matrix, const DenseVector& mean_vector,
//                    RandomAccessIterator begin, RandomAccessIterator end,
//                    FeatureVectorCallback callback, IndexType dimension)
//{
//    timed_context context("Data projection");
//    
//    DenseVector current_vector(dimension);
//    DenseVector current_vector_subtracted_mean(dimension);
//    
//    DenseMatrix embedding = DenseMatrix::Zero((end-begin),projection_matrix.cols());
//    
//    for (RandomAccessIterator iter=begin; iter!=end; ++iter)
//    {
//        callback.vector(*iter,current_vector);
//        current_vector_subtracted_mean = current_vector - mean_vector;
//        embedding.row(iter-begin) = projection_matrix.transpose()*current_vector_subtracted_mean;
//    }
//    
//    return embedding;
//}

//inline Matrix randomProjection(const Matrix &data, const int target_dimension) {
//    size_t current_dimension = data.cols();
//    
//    Matrix projection_matrix = gaussian_projection_matrix((int)current_dimension, target_dimension);
//    Matrix mean_vector = data.mean();// compute_mean(begin, end, features, current_dimension);
//    
//    tapkee::ProjectingFunction projecting_function(new tapkee::MatrixProjectionImplementation(projection_matrix,mean_vector));
//    
//    return TapkeeOutput(project(projection_matrix, mean_vector, begin, end, features, current_dimension), projecting_function);
//}

#endif
