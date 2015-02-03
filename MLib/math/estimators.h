//
//  estimators.h
//  TripSafety
//
//  Created by Iaroslav Omelianenko on 1/31/15.
//  Copyright (c) 2015 Danfoss. All rights reserved.
//

#ifndef TripSafety_estimators_h
#define TripSafety_estimators_h

#include "statistics.h"
#include "matrix.h"

namespace nologin {
    namespace math {
        /** The small deviation allowed in double comparisons. */
        const static double SMALL = 1e-6;
        
        class Estimator {
        public:
            /**
             * Get a probability estimate for a value
             *
             * @param data the value to estimate the probability of
             * @return the estimated probability of the supplied value
             */
            virtual double getProbability(double data) const = 0;
            
            /**
             * Add a new data value to the current estimator.
             *
             * @param data the new data value
             * @param weight the weight assigned to the data value
             */
            virtual void addValue(double data, double weight) = 0;
        };
        
        /**
         * Simple probability estimator that places a single Poisson distribution
         * over the observed values.
         */
        class PoissonEstimator : public Estimator {
            /** The number of values seen */
            double m_NumValues;
            /** The sum of the values seen */
            double m_SumOfValues;
            /**
             * The average number of times
             * an event occurs in an interval.
             */
            double m_Lambda;
        public:
            
            /**
             * Add a new data value to the current estimator.
             *
             * @param data the new data value
             * @param weight the weight assigned to the data value
             */
            void addValue(double data, double weight) {
                m_NumValues += weight;
                m_SumOfValues += data * weight;
                if (m_NumValues != 0) {
                    m_Lambda = m_SumOfValues / m_NumValues;
                }
            }
            
            /**
             * Get a probability estimate for a value
             *
             * @param data the value to estimate the probability of
             * @return the estimated probability of the supplied value
             */
            double getProbability(double data) const {
                double p = Poisson(data);
                return p;
            }
            
        private:
            /**
             * Calculates the log factorial of a number.
             *
             * @param x input number.
             * @return log factorial of x.
             */
            double logFac(double x) const {
                double result = 0;
                for (double i = 2; i <= x; i++) {
                    result += std::log(i);
                }
                return result;
            }
            
            /**
             * Returns value for Poisson distribution
             *
             * @param x the argument to the kernel function
             * @return the value for a Poisson kernel
             */
            double Poisson(double x) const {
                return std::exp(-m_Lambda + (x * std::log(m_Lambda)) - logFac(x));
            }
        };
        
        /**
         * Simple probability estimator that places a single normal distribution
         * over the observed values.
         *
         * @author Len Trigg (trigg@cs.waikato.ac.nz)
         * @version $Revision: 5540 $
         */
        class NormalEstimator : public Estimator {
            /** The sum of the weights */
            double m_SumOfWeights;
            /** The sum of the values seen */
            double m_SumOfValues;
            /** The sum of the values squared */
            double m_SumOfValuesSq;
            /** The current mean */
            double m_Mean;
            /** The current standard deviation */
            double m_StandardDev;
            /** The precision of numeric values ( = minimum std dev permitted) */
            double m_Precision;
            
        public:
            /**
             * Constructor that takes a precision argument.
             *
             * @param precision the precision to which numeric values are given. For
             * example, if the precision is stated to be 0.1, the values in the
             * interval (0.25,0.35] are all treated as 0.3.
             */
            NormalEstimator(double precision) {
                m_Precision = precision;
                // Allow at most 3 sd's within one interval
                m_StandardDev = m_Precision / (2 * 3);
            }
            
            /**
             * Add a new data value to the current estimator.
             *
             * @param data the new data value
             * @param weight the weight assigned to the data value
             */
            void addValue(double data, double weight) {
                if (weight == 0) {
                    return;
                }
                data = round(data);
                m_SumOfWeights += weight;
                m_SumOfValues += data * weight;
                m_SumOfValuesSq += data * data * weight;
                
                if (m_SumOfWeights > 0) {
                    m_Mean = m_SumOfValues / m_SumOfWeights;
                    double stdDev = std::sqrt(std::abs(m_SumOfValuesSq - m_Mean * m_SumOfValues) / m_SumOfWeights);
                    // If the stdDev ~= 0, we really have no idea of scale yet,
                    // so stick with the default. Otherwise...
                    if (stdDev > 1e-10) {
                        m_StandardDev = std::max(m_Precision / (2 * 3),
                                                 // allow at most 3sd's within one interval
                                                 stdDev);
                    }
                }
            }
            
            /**
             * Get a probability estimate for a value
             *
             * @param data the value to estimate the probability of
             * @return the estimated probability of the supplied value
             */
            double getProbability(double data) const {
                
                data = round(data);
                double zLower = (data - m_Mean - (m_Precision / 2)) / m_StandardDev;
                double zUpper = (data - m_Mean + (m_Precision / 2)) / m_StandardDev;
                
                double pLower = normalProbability(zLower);
                double pUpper = normalProbability(zUpper);
                
                double p = pUpper - pLower;
                return p;
            }
            
        private:
            /**
             * Round a data value using the defined precision for this estimator
             *
             * @param data the value to round
             * @return the rounded data value
             */
            double round(double data) const {
                
                return std::rint(data / m_Precision) * m_Precision;
            }
        };
        
        /**
         * Simple kernel density estimator. Uses one gaussian kernel per observed
         * data value.
         */
        class KernelEstimator : public Estimator {
            /** Vector containing all of the values seen */
            std::vector<double> m_Values;
            /** Vector containing the associated weights */
            std::vector<double> m_Weights;
            /** Number of values stored in m_Weights and m_Values so far */
            int m_NumValues;
            /** The sum of the weights so far */
            double m_SumOfWeights;
            /** The standard deviation */
            double m_StandardDev;
            /** The precision of data values */
            double m_Precision;
            /** Whether we can optimise the kernel summation */
            bool m_AllWeightsOne;
            /** Maximum percentage error permitted in probability calculations */
            double MAX_ERROR = 0.01;
            
        public:
            /**
             * Constructor that takes a precision argument.
             *
             * @param precision the  precision to which numeric values are given. For
             * example, if the precision is stated to be 0.1, the values in the
             * interval (0.25,0.35] are all treated as 0.3.
             */
            KernelEstimator(double precision) {
                
                m_Values.resize(50, 0);
                m_Weights.resize(50, 0);
                
                m_NumValues = 0;
                m_SumOfWeights = 0;
                m_AllWeightsOne = true;
                m_Precision = precision;
                // precision cannot be zero
                if (m_Precision < SMALL) m_Precision = SMALL;
                m_StandardDev = m_Precision / (2 * 3);
            }
            
            /**
             * Add a new data value to the current estimator.
             *
             * @param data the new data value
             * @param weight the weight assigned to the data value
             */
            void addValue(double data, double weight) {
                
                if (weight == 0) {
                    return;
                }
                data = round(data);
                int insertIndex = findNearestValue(data);
                if ((m_NumValues <= insertIndex) || (m_Values[insertIndex] != data)) {
                    if (m_NumValues < m_Values.size()) {
                        int left = m_NumValues - insertIndex;
                        std::copy(m_Values.begin() + insertIndex, m_Values.begin() + insertIndex + left, m_Values.begin() + insertIndex + 1);
                        std::copy(m_Weights.begin() + insertIndex, m_Weights.begin() + insertIndex + left, m_Weights.begin() + insertIndex + 1);
                        m_Values[insertIndex] = data;
                        m_Weights[insertIndex] = weight;
                        m_NumValues++;
                    } else {
                        std::vector<double>newValues(m_Values.size() * 2, 0);
                        std::vector<double>newWeights(m_Values.size() * 2, 0);
                        int left = m_NumValues - insertIndex;
                        std::copy(m_Values.begin(), m_Values.begin() + insertIndex, newValues.begin());
                        std::copy(m_Weights.begin(), m_Weights.begin() + insertIndex, newWeights.begin());
                        newValues[insertIndex] = data;
                        newWeights[insertIndex] = weight;
                        std::copy(m_Values.begin() + insertIndex, m_Values.begin() + insertIndex + left, newValues.begin() + insertIndex + 1);
                        std::copy(m_Weights.begin() + insertIndex, m_Weights.begin() + insertIndex + left, newWeights.begin() + insertIndex + 1);
                        m_NumValues++;
                        m_Values = newValues;
                        m_Weights = newWeights;
                    }
                    if (weight != 1) {
                        m_AllWeightsOne = false;
                    }
                } else {
                    m_Weights[insertIndex] += weight;
                    m_AllWeightsOne = false;
                }
                m_SumOfWeights += weight;
                double range = m_Values[m_NumValues - 1] - m_Values[0];
                if (range > 0) {
                    m_StandardDev = std::max(range / std::sqrt(m_SumOfWeights),
                                             // allow at most 3 sds within one interval
                                             m_Precision / (2 * 3));
                }
            }
            
            /**
             * Get a probability estimate for a value.
             *
             * @param data the value to estimate the probability of
             * @return the estimated probability of the supplied value
             */
            double getProbability(double data) const {
                
                double delta = 0, sum = 0, currentProb = 0;
                double zLower = 0, zUpper = 0;
                if (m_NumValues == 0) {
                    zLower = (data - (m_Precision / 2)) / m_StandardDev;
                    zUpper = (data + (m_Precision / 2)) / m_StandardDev;
                    return (normalProbability(zUpper) - normalProbability(zLower));
                }
                double weightSum = 0;
                int start = findNearestValue(data);
                for (int i = start; i < m_NumValues; i++) {
                    delta = m_Values[i] - data;
                    zLower = (delta - (m_Precision / 2)) / m_StandardDev;
                    zUpper = (delta + (m_Precision / 2)) / m_StandardDev;
                    currentProb = (normalProbability(zUpper) - normalProbability(zLower));
                    sum += currentProb * m_Weights[i];
                    weightSum += m_Weights[i];
                    if (currentProb * (m_SumOfWeights - weightSum) < sum * MAX_ERROR) {
                        break;
                    }
                }
                for (int i = start - 1; i >= 0; i--) {
                    delta = m_Values[i] - data;
                    zLower = (delta - (m_Precision / 2)) / m_StandardDev;
                    zUpper = (delta + (m_Precision / 2)) / m_StandardDev;
                    currentProb = (normalProbability(zUpper) - normalProbability(zLower));
                    sum += currentProb * m_Weights[i];
                    weightSum += m_Weights[i];
                    if (currentProb * (m_SumOfWeights - weightSum) < sum * MAX_ERROR) {
                        break;
                    }
                }
                double p = sum / m_SumOfWeights;
                return p;
            }
            
            
        private:
            /**
             * Execute a binary search to locate the nearest data value
             *
             * @param the data value to locate
             * @return the index of the nearest data value
             */
            int findNearestValue(double key) const {
                int low = 0;
                int high = m_NumValues;
                int middle = 0;
                while (low < high) {
                    middle = (low + high) / 2;
                    double current = m_Values[middle];
                    if (current == key) {
                        return middle;
                    }
                    if (current > key) {
                        high = middle;
                    } else if (current < key) {
                        low = middle + 1;
                    }
                }
                return low;
            }
            
            /**
             * Round a data value using the defined precision for this estimator
             *
             * @param data the value to round
             * @return the rounded data value
             */
            double round(double data) {
                return std::rint(data / m_Precision) * m_Precision;
            }
        };
        
        /**
         * Simple symbolic probability estimator based on symbol counts.
         */
        class DiscreteEstimator : public Estimator {
            /** Hold the counts */
            std::vector<double> m_Counts;
            /** Hold the sum of counts */
            double m_SumOfCounts;
            
        public:
            /**
             * Constructor
             *
             * @param numSymbols the number of possible symbols (remember to include 0)
             * @param laplace if true, counts will be initialised to 1
             */
            DiscreteEstimator(int numSymbols, bool laplace) {
                if (laplace) {
                    m_Counts.resize(numSymbols, 1);
                    m_SumOfCounts = (double)numSymbols;
                } else {
                    m_Counts.resize(numSymbols, 0);
                    m_SumOfCounts = 0;
                }
            }
            
            /**
             * Constructor
             *
             * @param nSymbols the number of possible symbols (remember to include 0)
             * @param fPrior value with which counts will be initialised
             */
            DiscreteEstimator(int nSymbols, double fPrior) {
                m_Counts.resize(nSymbols, 0);
                for(int iSymbol = 0; iSymbol < nSymbols; iSymbol++) {
                    m_Counts[iSymbol] = fPrior;
                }
                m_SumOfCounts = fPrior * (double) nSymbols;
            }
            
            /**
             * Add a new data value to the current estimator.
             *
             * @param data the new data value
             * @param weight the weight assigned to the data value
             */
            void addValue(double data, double weight) {
                if ((int)data < m_Counts.size()) {
                    m_Counts[(int)data] += weight;
                    m_SumOfCounts += weight;
                }
            }
            
            /**
             * Get a probability estimate for a value
             *
             * @param data the value to estimate the probability of
             * @return the estimated probability of the supplied value
             */
            double getProbability(double data) const {
                if (m_SumOfCounts == 0) {
                    return 0;
                }
                double p = (double)m_Counts[(int)data] / m_SumOfCounts;
                return p;
            }
        };
        
        class MahalanobisEstimator : public Estimator{
            /** The inverse of the covariance matrix */
            Matrix *m_CovarianceInverse;
            /** The determinant of the covariance matrix */
            double m_Determinant;
            /** The difference between the conditioning value and the conditioning mean */
            double m_ConstDelta;
            /** The mean of the values */
            double m_ValueMean;
            
        public:
            MahalanobisEstimator(const Matrix &covariance, const double constDelta, const double valueMean) {
                Assert((covariance.rows() == 2) && (covariance.cols() == 2),
                       "Wrong covariance matrix dimensions! Rows: %lu, cols: %lu", covariance.rows(), covariance.cols());
                m_CovarianceInverse = NULL;
                
                double a = covariance(0, 0);
                double b = covariance(0, 1);
                double c = covariance(1, 0);
                double d = covariance(1, 1);
                if (a == 0) {
                    a = c; c = 0;
                    double temp = b;
                    b = d; d = temp;
                }
                if (a == 0) {
                    return;
                }
                double denom = d - c * b / a;
                if (denom == 0) {
                    return;
                }
                
                m_Determinant = covariance(0, 0) * covariance(1, 1) - covariance(1, 0) * covariance(0, 1);
                m_CovarianceInverse = new Matrix(2, 2);
                (*m_CovarianceInverse)(0, 0) = 1.0 / a + b * c / a / a / denom;
                (*m_CovarianceInverse)(0, 1) = -b / a / denom;
                (*m_CovarianceInverse)(1, 0) = -c / a / denom;
                (*m_CovarianceInverse)(1, 1) = 1.0 / denom;
                m_ConstDelta = constDelta;
                m_ValueMean = valueMean;
                
            }
            
            void addValue(double data, double weight) {
                
            }
            
            double getProbability(double data) const {
                
                double delta = data - m_ValueMean;
                if (m_CovarianceInverse == NULL) {
                    return 0;
                }
                return normalKernel(delta);
            }
            
        private:
            double normalKernel(double x) const {
                Matrix thisPoint(1, 2);
                thisPoint(0, 0) = x;
                thisPoint(0, 1) = m_ConstDelta;
                
                thisPoint *= (*m_CovarianceInverse);
                thisPoint *= thisPoint.transpose();
                double val = exp(thisPoint(0, 0) / 2) / (sqrt(M_PI * 2) * m_Determinant);
                
                //        exp(-thisPoint.times(m_CovarianceInverse).times(thisPoint.transpose()).get(0, 0)/ 2) / (sqrt(M_PI * 2) * m_Determinant);
                return val;
            }
            
        };
        class NNConditionalEstimator {
            /** Vector containing all of the values seen */
            VD m_Values;
            /** Vector containing all of the conditioning values seen */
            VD m_CondValues;
            /** Vector containing the associated weights */
            VD m_Weights;
            /** The sum of the weights so far */
            double m_SumOfWeights;
            /** Current Conditional mean */
            double m_CondMean;
            /** Current Values mean */
            double m_ValueMean;
            /** Current covariance matrix */
            Matrix *m_Covariance;
            /** Whether we can optimise the kernel summation */
            bool m_AllWeightsOne = true;
            
        public:
            void addValue(double data, double given, double weight) {
                
                size_t insertIndex = findNearestPair(given, data);
                
                
                if ((m_Values.size() <= insertIndex) || (m_CondValues[insertIndex] != given) || (m_Values[insertIndex] != data)) {
                    m_CondValues.insert(m_CondValues.begin() + insertIndex, given);
                    m_Values.insert(m_Values.begin() + insertIndex, data);
                    m_Weights.insert(m_Weights.begin() + insertIndex, weight);
                    if (weight != 1) {
                        m_AllWeightsOne = false;
                    }
                } else {
                    double newWeight = m_Weights[insertIndex];
                    newWeight += weight;
                    m_Weights[insertIndex] = newWeight;
                    m_AllWeightsOne = false;
                }
                m_SumOfWeights += weight;
                
                // Invalidate any previously calculated covariance matrix
                delete m_Covariance;
                m_Covariance = NULL;
            }
            
            double getProbability(double data, double given) {
                if (m_Covariance == NULL) {
                    calculateCovariance();
                }
                MahalanobisEstimator estimator(*m_Covariance, given - m_CondMean, m_ValueMean);
                return estimator.getProbability(data);
            }
            
        private:
            size_t findNearestPair(double key, double secondaryKey) {
                size_t low = 0;
                size_t high = m_CondValues.size();
                size_t middle = 0;
                while (low < high) {
                    middle = (low + high) / 2;
                    double current = m_CondValues[middle];
                    if (current == key) {
                        double secondary = m_Values[middle];
                        if (secondary == secondaryKey) {
                            return middle;
                        }
                        if (secondary > secondaryKey) {
                            high = middle;
                        } else if (secondary < secondaryKey) {
                            low = middle + 1;
                        }
                    }
                    if (current > key) {
                        high = middle;
                    } else if (current < key) {
                        low = middle + 1;
                    }
                }
                return low;
            }
            
            void calculateCovariance() {
                double sumValues = 0, sumConds = 0;
                for(int i = 0; i < m_Values.size(); i++) {
                    sumValues += m_Values[i] * m_Weights[i];
                    sumConds += m_CondValues[i] * m_Weights[i];
                }
                m_ValueMean = sumValues / m_SumOfWeights;
                m_CondMean = sumConds / m_SumOfWeights;
                double c00 = 0, c01 = 0, c10 = 0, c11 = 0;
                for(int i = 0; i < m_Values.size(); i++) {
                    double x = m_Values[i];
                    double y = m_CondValues[i];
                    double weight = m_Weights[i];
                    
                    Printf("x: %f, y: %f, weight: %f\n", x, y, weight);
                    c00 += (x - m_ValueMean) * (x - m_ValueMean) * weight;
                    c01 += (x - m_ValueMean) * (y - m_CondMean) * weight;
                    c11 += (y - m_CondMean) * (y - m_CondMean) * weight;
                }
                c00 /= (m_SumOfWeights - 1.0);
                c01 /= (m_SumOfWeights - 1.0);
                c10 = c01;
                c11 /= (m_SumOfWeights - 1.0);
                m_Covariance = new Matrix(2, 2);
                (*m_Covariance)(0, 0) = c00;
                (*m_Covariance)(0, 1) = c01;
                (*m_Covariance)(1, 0) = c10;
                (*m_Covariance)(1, 1) = c11;
            }
        };
    }
}
#endif
