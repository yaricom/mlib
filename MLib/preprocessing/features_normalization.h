//
//  features_normalization.h
//  MLib
//
//  Created by Iaroslav Omelianenko on 1/28/15.
//  Copyright (c) 2015 yaric. All rights reserved.
//

#ifndef MLib_features_normalization_h
#define MLib_features_normalization_h

#include "matrix.h"
#include "random.h"

namespace nologin {
    namespace preprocessing {
        
        /**
         * Scales this matrix to have all samles scaled to fit range [min, max] per features vectors
         *
         * X_std = (X - X.min) / (X.max - X.min)
         * X_scaled = X_std * (max - min) + min
         *
         * @param min the minimal range limit inclusive
         * @param max the maximal range limit inclusive
         * @param indices the column's indices for processing
         * @param mat the matrix to be processed
         */
        void scaleMinMax(const double min, const double max, const VI &indices, Matrix &mat) {
            size_t n = mat.cols();
            size_t m = mat.rows();
            
            size_t ind_size = indices.size();
            assert(ind_size <= n);
            
            // fin min/max per sample per feature
            VD fMins(n, numeric_limits<double>().max()), fMaxs(n, numeric_limits<double>().min());
            for (int row = 0 ; row < m; row++) {
                for (int col : indices) {
                    assert(col < n);
                    fMaxs[col] = std::max<double>(fMaxs[col], mat[row][col]);
                    fMins[col] = std::min<double>(fMins[col], mat[row][col]);
                }
            }
            
            // find X scaled
            for (int row = 0 ; row < m; row++) {
                for (int col : indices) {
                    double X = mat[row][col], X_min = fMins[col], X_max = fMaxs[col];
                    double X_std = (X - X_min) / (X_max - X_min);
                    mat[row][col] = X_std * (max - min) + min;
                }
            }
        }
        
        /**
         * Scales this matrix to have all samles scaled to fit range [min, max] per features vectors
         *
         * @param min the minimal range limit inclusive
         * @param max the maximal range limit inclusive
         * @param mat the matrix to be processed
         */
        void scaleMinMax(const double min, const double max, Matrix &mat) {
            size_t n = mat.cols();
            // create indices array
            VI indices(n, -1);
            for (int j = 0; j < n; j++) {
                indices[j] = j;
            }
            
            scaleMinMax(min, max, indices, mat);
        }
        
        /**
         * Scales this matrix to have all values centered arround zero with standard deviation = 1
         *
         * @param indices the column's indices for processing
         * @param mat the matrix to be processed
         */
        void stdScale(const VI &indices, Matrix &mat) {
            size_t n = mat.cols();
            size_t m = mat.rows();
            
            size_t ind_size = indices.size();
            assert(ind_size <= n);
            
            Vector meanV = mat.mean();
            Vector stdV = mat.stdev(1);
            
            for (int row = 0 ; row < m; row++) {
                for (int col : indices) {
                    double X = mat[row][col];
                    mat[row][col] = (X - meanV[col]) / stdV[col];
                }
            }
        }
        
        /**
         * Scales this matrix to have all values centered arround zero with standard deviation = 1
         * @param mat the matrix to be processed
         */
        void stdScale(Matrix &mat) {
            size_t n = mat.cols();
            
            // create indices array
            VI indices(n, -1);
            for (int j = 0; j < n; j++) {
                indices[j] = j;
            }
            stdScale(indices, mat);
        }
        
        /**
         * Correct outliers in the specified columns.
         * @param indices The columns indices.
         * @param of The factor for outlier detection.
         * @param evf The factor for extreme values detection.
         * @param mat the matrix to be processed
         */
        void correctOutliers(const VI &indices, Matrix &mat, const double of = 1.5, const double evf = 3) {
            size_t indCount = indices.size();
            assert(indCount < mat.cols());
            
            int	half, quarter;
            double q1, q2, q3;
            
            VD upperExtremeValue(indCount, 0);
            VD upperOutlier(indCount, 0);
            VD lowerOutlier(indCount, 0);
            VD lowerExtremeValue(indCount, 0);
            VD median(indCount, 0);
            VD IQR(indCount, 0);
            
            //
            // find inetrquartile distances
            //
            VD values(mat.rows(), 0);
            for (int j = 0; j < indCount; j++) {
                mat.columnToArray(indices[j], values);
                std::sort(values.begin(), values.end());
                
                // determine indices
                half = (int)values.size() / 2;
                quarter = half / 2;
                
                // find quartiles
                if (values.size() % 2 == 1) {
                    // odd
                    q2 = values[half];
                } else {
                    // even
                    q2 = (values[half] + values[half - 1]) / 2;
                }
                
                if (half % 2 == 1) {
                    // odd
                    q1 = values[quarter];
                    q3 = values[values.size() - quarter - 1];
                }
                else {
                    // even
                    q1 = (values[quarter] + values[quarter - 1]) / 2;
                    q3 = (values[values.size() - quarter - 1] + values[values.size() - quarter]) / 2;
                }
                
                // determine thresholds and other values
                median[j] = q2;
                IQR[j] = q3 - q1;
                upperExtremeValue[j] = q3 + evf * IQR[j];
                upperOutlier[j] = q3 + of * IQR[j];
                lowerOutlier[j] = q1 - of * IQR[j];
                lowerExtremeValue[j] = q1 - evf * IQR[j];
            }
            
            //
            // Correct input matrix
            //
            // obtain a seed from the system clock:
            unsigned int seed = (unsigned int)std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine generator(seed);
            std::uniform_real_distribution<double> rnd(0.0, 1.0);
            size_t m = mat.rows();
            for (int row = 0 ; row < m; row++) {
                for (int j = 0; j < indCount; j++) {
                    double val = mat[row][indices[j]];
                    if (val < lowerExtremeValue[j]) {
                        mat[row][indices[j]] = lowerOutlier[j];
                    } else if (val < lowerOutlier[j]) {
                        double lowerRange = median[j] - lowerOutlier[j];
                        mat[row][indices[j]] = median[j] - rnd(generator) * lowerRange;
                    } else if (val > upperExtremeValue[j]) {
                        mat[row][indices[j]] = upperOutlier[j];
                    } else if (val > upperOutlier[j]){
                        double upperRange = upperOutlier[j] - median[j];
                        mat[row][indices[j]] = median[j] + rnd(generator) * upperRange;
                    }
                }
            }
        }
        
        VD expandRow(const VD &row) {
            size_t N = row.size();
            
            VD rv, rvM, rvD;
            for (int i = 0; i < N; i++) {
                rv.push_back(row[i]);
                
                for (int j = 0; j < i; j++) {
                    // find multiplication
                    rvM.push_back(row[i] * row[j]);
                    // find division
                    double v = row[j] ? row[i] / row[j] : 0;
                    rvD.push_back(v);
                }
            }
            
            // append vectors
            rv.reserve(rv.size() + rvM.size() + rvD.size());
            rv.insert(rv.end(), rvM.begin(), rvM.end());
            rv.insert(rv.end(), rvD.begin(), rvD.end());
            
            return rv;
        }
        
        Matrix& expandMatrixLineary(const Matrix &mat) {
            VVD rv;
            for (int i = 0; i < mat.rows(); i++) {
                rv.push_back(expandRow(mat[i]));
            }
            
            Matrix *ret = new Matrix(rv);
            return *ret;
        }
        
    }
}

#endif
