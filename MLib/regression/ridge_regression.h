//
//  ridge_regression.h
//  TripSafety
//
//  Created by Iaroslav Omelianenko on 2/10/15.
//  Copyright (c) 2015 nologin. All rights reserved.
//

#ifndef TripSafety_ridge_regression_h
#define TripSafety_ridge_regression_h

/*
 The data is read from tab-separated ASCI file rg.txt; each column is
 saved into auxilary files rg_xx.txt.
 
 File rg.txt: data set; column 0 = dependent var
 Files rg_xx.txt are copies of individual columns
 File rg_err.txt: residual error per observation
 File rg_log.txt: program logs: this file contains all results (regression
 coefficients, empirical distribution for regression coefficients, etc.)
 
 The structure of rg_log.txt easily allows for parsing and further
 statistical processing:
 
 col 1: xxxxx_yy where xxxxx is the name of the procedure being run
 and yy is the iteration.
 col 2: iteration in the Regression call; set VERBOSE to -1 if you
 don't want this level of detail - then only results obtained
 at the final iteration are saved in rg_log.txt
 col 3: variable optimized at current iteration
 col 4: lambda (internal variable)
 col 5: reduction in standard deviation of error achieved at current
 iteration; this value is between 0 and 1; 0 corresponds to all
 regression coefficients set to zero; 1 means perfect fit
 cols 6, 7, 8 etc.: current value of regression coefficients
 
 *param: regression coefficients
 *param_seed: initial regression coefficients (usually 0)
 obs: number of observations
 var: number of columns in input file
 MSE: mean squared error
 MSE_init: MSE when regression coefficients are 0
 init: flag to indicate that Regress_init has been called.
 seed: flag to indicate that Init_param has been called. */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "utils.h"
#include "defines.h"

using namespace std;

namespace nologin {
    namespace regression {
        
        using namespace utils;
        
#define VERBOSE 0  /* 0 for detailed logs, -1 for summary */
        
        struct RRConfig {
            bool useBootstrap;
            int bootstrapSamples;
            int bootstrapMode;
            int bootstrapIterations;
            int regressionIterations;
            int regressionMode;
            
            int approxSamples;
            int approxMode;
            int approxIterations;
        };
        
        class RidgeRegression {
            // regression coefficients
            VD regrCoef;
            // initial regression coefficients (usually 0)
            VD regrCoefSeed;
            // residual errors per observation
            VD obsErrors;
            
            // the OOB errors calculated if bootstrap used
            VD oobErrors;
            
            // the configuration
            RRConfig config;
            
            // the number of observations/samples
            size_t obs = -1;
            // the number of features per observation
            size_t var = -1;
            
            // flag to indicate that data structures initialize
            int init = -1;
            // flag to indicate that approximate regression coeff was calculated
            int seed = -1;
            
            // the calculated MSE after regression complete
            double MSE = 0;
            // the initial MSE before regresion starts
            double MSE_init = 0;
            
        public:
            RidgeRegression(RRConfig config) : config(config) {}
            
            /**
             * Start training regression tree by finding regression coefficients
             * @param train the vector of samples with features per column
             * @param check the vector of dependent variables per sample
             */
            void train(const VVD &train, const VD &check) {
                Assert(train.size() == check.size(), "Samples size should be equal to observations size");
                // find number of variables and observations
                var = train[0].size();
                obs = train.size();
                
                initDataStructures();
                if (config.useBootstrap) {
                    bootstrap(train, check, config.bootstrapSamples, config.bootstrapMode, config.bootstrapIterations);
                    regress(train, check, config.regressionMode, config.regressionIterations, 1);
                } else {
                    initFeaturesSeed(train, check, config.approxMode, config.approxSamples, config.approxIterations);
                    regress(train, check, config.regressionMode, config.regressionIterations, 1);
                }
            }
            
            /**
             * Do predict based on provided features
             * @param features the features of data sample to find DV
             * @return predicted dependent variable
             */
            double predict(const VD &features) {
                Assert(features.size() == var, "Features size should be equal to train features size, but was: %li", features.size());
                double res = 0;
                for (int i = 0; i < features.size(); i++) {
                    res += features[i] * regrCoef[i];
                }
                return res;
            }
            
        private:
            
            /**
             * Performs nvalid approximate regression with different starting
             * points to identify an initial approximate solution. This
             * solution will be used as starting point (seed) for the
             * regression.
             *
             * @param mode should be 2 (recommended) or 3; see Regression for details
             * @param nvalid: number of regressions to perform
             * @param niter: number of iterations to use in each regression; should be small here (<10)
             * @param mode: should be 2 (recommended) or 3; see Regression for description
             * Output (global): seed is set to 1 to indicate that Regression must use regrCoefSeed
             * Output: regrCoefSeed: regression coefficients to be used as starting point in Regression
             */
            void initFeaturesSeed(const VVD &train, const VD &check, const long mode, const long nvalid, const long niter) {
                // initialize data structure
                if (seed==1) {
                    regrCoefSeed.clear();
                } else {
                    seed=1;
                }
                regrCoefSeed.resize(var, 0);
                
                double nsvar, nsvar_max = -1;
                for (int n = 0; n < nvalid; n++) {
                    nsvar = regress(train, check, mode, niter, 0);
                    if (nsvar > nsvar_max) {
                        nsvar_max = nsvar;
                        for (int k = 0; k < var; k++) {
                            regrCoefSeed[k] = regrCoef[k];
                        }
                    }
                }
            }
            
            /**
             * Compute empirical distribution for estimated regression coefficients
             * and reduction in standard deviation of error. Could also be used to
             * compute empirical distribution of error for each observation, to
             * detect outliers.
             *
             * Input (global): var, obs, init (must be set to 1 by using initDataStructures() first)
             * @param niter: number of iterations to use to compute regression coefficients; see Regression.
             * @param nsample: number of samples; Boostrap performs one regression on each sample
             * @param mode: see Regression for details.
             *
             * Output (global): regrCoef: the optimal found regression coefficients
             */
            void bootstrap(const VVD &train, const VD &check, const long nsample, const long mode, const long niter) {
                /* need to run Regress_init first if the data set is new */
                Assert(init == 1, "Must run initDataStructures() fisrt.\n");
                
                // initialize data structure
                if (seed==1) {
                    regrCoefSeed.clear();
                } else {
                    seed=1;
                }
                regrCoefSeed.resize(var, 0);
                
                VVD bootTrain, bootTest;
                VD bootCheck, bootCheckTest;
                
                // the best found regression coefficients
                VD bestCoef(var, 0);
                
                // array to store selected indices
                VI pick(obs, 0);
                int k, idx, m, oi;
                double oobMSE, oobME, minOOBMSE = numeric_limits<double>::max();
                for (int n = 0; n < nsample; n++) {
                    // pick up random observations
                    for (k = 0; k < obs; k++) { pick[k] = 0; }
                    for (k = 0; k < obs; k++) {
                        idx = rand() % obs;
                        pick[idx]++;
                    }
                    
                    // create subsample
                    for (m = 0; m < obs; m++) {
                        if (pick[m] > 0) {
                            // save pick[m] copies of row in data
                            for (k = 0; k < pick[m]; k++) {
                                bootTrain.push_back(train[m]);
                                bootCheck.push_back(check[m]);
                            }
                        } else {
                            // save current row as test sample
                            bootTest.push_back(train[m]);
                            bootCheckTest.push_back(check[m]);
                        }
                    }
                    // do regression
                    initDataStructures();
                    regress(bootTrain, bootCheck, mode, niter, 0);
                    
                    // find OOB error
                    oobMSE = 0;
                    for (oi = 0; oi < bootTest.size(); oi++) {
                        oobME = predict(bootTest[oi]) - bootCheckTest[oi];
                        oobMSE += oobME * oobME;
                    }
                    oobMSE = sqrt(oobMSE /(double)bootTest.size());
                    oobErrors.push_back(oobMSE);
                    if (oobMSE < minOOBMSE) {
                        oobMSE = minOOBMSE;
                        // store current regression coefficients
                        bestCoef.swap(regrCoef);
                    }
                }
                
                // store best regression coefficients as train results
                //        regrCoef.swap(bestCoef);
                regrCoefSeed.swap(bestCoef);
            }
            
            /**
             * Performs regression. Must run initDataStructures() first to compute obs,
             * var, initialize init and allocate memory to param. The Regress function
             * returns the reduction in standard deviation of error; this value is
             * between 0 and 1; 0 corresponds to all regression coefficients set to
             * zero; 1 means perfect fit. Regress does not compute confidence
             * intervals for the coefficients (use Boostrap for this purpose).
             *
             * Input (global): var, obs, seed, init
             * @param seed_flag if 1 use initial regression coefficients param_seed
             * to start iterative regression procedure; if seed_flag = 1 then
             * Param_init must be called first to intialize  param_seed and seed.
             *
             * @param niter number of iterations to use to compute regression
             * coefficients; if convergence is slow or erratic then increase
             * niter and perform validation tests with Validation.
             *
             * @param mode determines the type of algorithm used for regression:
             *  mode = 0: visits each variable sequentially starting with first
             *  variable; useful when variables are pre-sorted in such a
             *  way that the first few variables explain most of the
             *  variance.
             *  mode = 1: visits variables in random order; should be the default mode
             *  mode = 2: same as mode = 1, but lambda not randomized
             *
             * @return reduction in standard deviation of error.
             * Output (global): regrCoef: regression coefficients
             */
            double regress(const VVD &features, const VD &dv, const long mode, const long niter, const long seed_flag) {
                /* need to run Regress_init first if the data set is new */
                Assert(init == 1, "Must run initDataStructures() fisrt.\n");
                
                int k, l, i, iter;
                double xd, sp, lambda, val, e_new, resvar = 0;
                
                // calculate initial MSE
                MSE_init = 0;
                for (k = 0; k < obs; k++) {
                    val = dv[k];
                    MSE_init += val * val;
                    obsErrors[k] = val;
                }
                MSE_init = sqrt(MSE_init / obs);
                
                // clear regression coefficients
                for (k = 1; k < var; k++) {
                    regrCoef[k] = 0;
                }
                
                // if seed=1 uses initial regressors
                if (seed_flag == 1) {
                    Assert(seed == 1, "Must run initFeaturesSeed() first.\n");
                    
                    for (k = 0; k < var; k++) {
                        regrCoef[k] = regrCoefSeed[k];
                    }
                    
                    for (i = 0; i < obs; i++) {
                        e_new = 0;
                        for (k = 0; k < var; k++) {
                            val = features[i][k];
                            e_new += regrCoef[k] * val;
                        }
                        // find error
                        obsErrors[i] = dv[i] - e_new;
                    }
                }
                
                /*
                 regression
                 */
                for (iter = 0; iter < niter; iter++) {
                    if (mode == 0 || mode == 3) {
                        // 0: visits each variable sequentially starting with first
                        l = iter % var;
                    } else {
                        // 1: visits variables in random order; should be the default mode
                        // 2: same as mode = 1, but lambda not randomized
                        l = rand() % var;
                    }
                    
                    xd = 0; sp = 0;
                    for (k = 0; k < obs; k++) {
                        xd += features[k][l] * features[k][l];
                        sp += features[k][l] * obsErrors[k];
                    }
                    Assert(xd != 0, "Empty column found at index: %i", l);
                    
                    lambda = sp / xd;
                    if (mode == 1) {
                        lambda = lambda * rand() / (double)RAND_MAX;
                    }
                    
                    // update error
                    MSE = 0;
                    for (k = 0; k < obs; k++) {
                        e_new = obsErrors[k] - lambda * features[k][l];
                        MSE += e_new * e_new;
                        obsErrors[k] = e_new;
                    }
                    regrCoef[l] += lambda;
                    
                    /*
                     save results, compute resvar
                     */
                    
                    if (iter % 10 == VERBOSE || iter == niter-1) {
                        MSE = sqrt(MSE / obs);
                        resvar = 1 - MSE / MSE_init;
                        if (LOG_DEBUG) {
                            Printf("REGRESS %d\t%d\t%f\t%f\t", iter, l, lambda, resvar);
                            print(regrCoef);
                        }
                    }
                }
                return resvar;
            }
            
            /**
             * Reinitialize internal data structures.
             */
            void initDataStructures() {
                // clear data if needed
                if (init == 1) {
                    regrCoef.clear();
                    obsErrors.clear();
                } else  {
                    init = 1;
                }
                
                // resize internal data structures
                regrCoef.resize(var, 0);
                obsErrors.resize(obs, 0);
            }
            
        };
    }
}
#endif
