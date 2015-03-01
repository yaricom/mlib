//
//  gradient_boosting_tree.h
//  MLib
//
//  Created by Iaroslav Omelianenko on 1/29/15.
//  Copyright (c) 2015 yaric. All rights reserved.
//

#ifndef MLib_gradient_boosting_tree_h
#define MLib_gradient_boosting_tree_h

#include "matrix.h"
#include "random.h"

namespace nologin {
    namespace tree {
        
        using namespace math;
        using namespace utils;
        
        struct GBTConfig {
            double sampling_size_ratio = 0.5;
            double learning_rate = 0.01;
            int tree_number = 130;
            int tree_min_nodes = 10;
            int tree_depth = 3;
        };
        
        struct Node {
            double m_node_value;
            int m_feature_index;
            double m_terminal_left;
            double m_terminal_right;
            
            // Each non-leaf node has a left child and a right child.
            Node *m_left_child = NULL;
            Node *m_right_child = NULL;
            
            // Construction function
            Node(double value, int feature_index, double value_left, double value_right) :
            m_node_value(value), m_feature_index(feature_index), m_terminal_left(value_left), m_terminal_right(value_right) {}
            
        private:
            
            Node(Node const &); // non construction-copyable
            Node& operator=(Node const &); // non copyable
        };
        
        struct BestSplit {
            // the index of split feature
            int m_feature_index;
            // the calculated node value
            double m_node_value;
            // if false - split failed
            bool m_status;
            
            // construction function
            BestSplit() : m_feature_index(0.0), m_node_value(0.0), m_status(false) {}
        };
        
        struct SplitRes {
            VC<VD> m_feature_left;
            VC<VD> m_feature_right;
            double m_left_value = 0.0;
            double m_right_value = 0.0;
            VD m_obs_left;
            VD m_obs_right;
            
            // construction function
            SplitRes() : m_left_value(0.0), m_right_value(0.0) {}
        };
        
        struct ListData {
            double m_x;
            double m_y;
            
            ListData(double x, double y) : m_x(x), m_y(y) {}
            
            bool operator < (const ListData& str) const {
                return (m_x < str.m_x);
            }
        };
        
        typedef enum _TerminalType {
            AVERAGE, MAXIMAL
        }TerminalType;
        
        class RegressionTree {
        private:
            // class members
            int m_min_nodes;
            int m_max_depth;
            int m_current_depth;
            TerminalType m_type;
            
        public:
            // The root node
            Node *m_root = NULL;
            // The features importance per index
            VI features_importance;
            
            // construction function
            RegressionTree() : m_min_nodes(10), m_max_depth(3), m_current_depth(0), m_type(AVERAGE) {}
            
            // set parameters
            void setMinNodes(int min_nodes) {
                Assert(min_nodes > 3, "The number of terminal nodes is too small: %i", min_nodes);
                m_min_nodes = min_nodes;
            }
            
            void setDepth(int depth) {
                Assert(depth > 0, "Tree depth must be positive: %i", depth);
                m_max_depth = depth;
            }
            
            // get fit value
            double predict(const VD &feature_x) const{
                double re_res = 0.0;
                
                if (!m_root) {
                    // failed in building the tree
                    return re_res;
                }
                
                Node *current = m_root;
                
                while (true) {
                    // current node information
                    int c_feature_index = current->m_feature_index;
                    double c_node_value = current->m_node_value;
                    double c_node_left_value = current->m_terminal_left;
                    double c_node_right_value = current->m_terminal_right;
                    
                    if (feature_x[c_feature_index] < c_node_value) {
                        // we should consider left child
                        current = current->m_left_child;
                        
                        if (!current) {
                            re_res = c_node_left_value;
                            break;
                        }
                    } else {
                        // we should consider right child
                        current = current->m_right_child;
                        
                        if (!current) {
                            re_res = c_node_right_value;
                            break;
                        }
                    }
                }
                
                return re_res;
            }
            
            /*
             *  The method to build regression tree
             */
            void buildRegressionTree(const VC<VD> &samples_x, const VD &obs_y) {
                size_t samples_num = samples_x.size();
                
                Assert(samples_num == obs_y.size() && samples_num != 0,
                       "The number of samles does not match with the number of observations or the samples number is 0. Samples: %i", samples_num);
                
                Assert (m_min_nodes * 2 <= samples_num, "The number of samples is too small");
                
                size_t feature_dim = samples_x[0].size();
                features_importance.resize(feature_dim, 0);
                
                // build the regression tree
                buildTree(samples_x, obs_y);
            }
            
        private:
            
            /*
             *  The following function gets the best split given the data
             */
            BestSplit findOptimalSplit(const VC<VD> &samples_x, const VD &obs_y) {
                
                BestSplit split_point;
                
                if (m_current_depth > m_max_depth) {
                    return split_point;
                }
                
                size_t samples_num = samples_x.size();
                
                if (m_min_nodes * 2 > samples_num) {
                    // the number of observations in terminals is too small
                    return split_point;
                }
                size_t feature_dim = samples_x[0].size();
                
                
                double min_err = 0;
                int split_index = -1;
                double node_value = 0.0;
                
                // begin to get the best split information
                for (int loop_i = 0; loop_i < feature_dim; loop_i++){
                    // get the optimal split for the loop_index feature
                    
                    // get data sorted by the loop_i-th feature
                    VC<ListData> list_feature;
                    for (int loop_j = 0; loop_j < samples_num; loop_j++) {
                        list_feature.push_back(ListData(samples_x[loop_j][loop_i], obs_y[loop_j]));
                    }
                    
                    // sort the list
                    sort(list_feature.begin(), list_feature.end());
                    
                    // begin to split
                    double sum_left = 0.0;
                    double mean_left = 0.0;
                    int count_left = 0;
                    double sum_right = 0.0;
                    double mean_right = 0.0;
                    int count_right = 0;
                    double current_node_value = 0;
                    double current_err = 0.0;
                    
                    // initialize left
                    for (int loop_j = 0; loop_j < m_min_nodes; loop_j++) {
                        ListData fetched_data = list_feature[loop_j];
                        sum_left += fetched_data.m_y;
                        count_left++;
                    }
                    mean_left = sum_left / count_left;
                    
                    // initialize right
                    for (int loop_j = m_min_nodes; loop_j < samples_num; loop_j++) {
                        ListData fetched_data = list_feature[loop_j];
                        sum_right += fetched_data.m_y;
                        count_right++;
                    }
                    mean_right = sum_right / count_right;
                    
                    // calculate the current error
                    // err = ||x_l - mean(x_l)||_2^2 + ||x_r - mean(x_r)||_2^2
                    // = ||x||_2^2 - left_count * mean(x_l)^2 - right_count * mean(x_r)^2
                    // = constant - left_count * mean(x_l)^2 - right_count * mean(x_r)^2
                    // Thus, we only need to check "- left_count * mean(x_l)^2 - right_count * mean(x_r)^2"
                    current_err = -1 * count_left * mean_left * mean_left - count_right * mean_right * mean_right;
                    
                    // current node value
                    current_node_value = (list_feature[m_min_nodes].m_x + list_feature[m_min_nodes - 1].m_x) / 2;
                    
                    if (current_err < min_err && current_node_value != list_feature[m_min_nodes - 1].m_x) {
                        split_index = loop_i;
                        node_value = current_node_value;
                        min_err = current_err;
                    }
                    
                    // begin to find the best split point for the feature
                    for (int loop_j = m_min_nodes; loop_j <= samples_num - m_min_nodes - 1; loop_j++) {
                        ListData fetched_data = list_feature[loop_j];
                        double y = fetched_data.m_y;
                        sum_left += y;
                        count_left++;
                        mean_left = sum_left / count_left;
                        
                        
                        sum_right -= y;
                        count_right--;
                        mean_right = sum_right / count_right;
                        
                        
                        current_err = -1 * count_left * mean_left * mean_left - count_right * mean_right * mean_right;
                        // current node value
                        current_node_value = (list_feature[loop_j + 1].m_x + fetched_data.m_x) / 2;
                        
                        if (current_err < min_err && current_node_value != fetched_data.m_x) {
                            split_index = loop_i;
                            node_value = current_node_value;
                            min_err = current_err;
                        }
                        
                    }
                }
                // set the optimal split point
                if (split_index == -1) {
                    // failed to split data
                    return split_point;
                }
                split_point.m_feature_index = split_index;
                split_point.m_node_value = node_value;
                split_point.m_status = true;
                
                return split_point;
            }
            
            /*
             *  Split data into the left node and the right node based on the best splitting
             *  point.
             */
            SplitRes splitData(const VC<VD> &samples_x, const VD &obs_y, const BestSplit &best_split) {
                
                SplitRes split_res;
                
                int feature_index = best_split.m_feature_index;
                double node_value = best_split.m_node_value;
                
                size_t samples_count = obs_y.size();
                for (int loop_i = 0; loop_i < samples_count; loop_i++) {
                    VD ith_sample = samples_x[loop_i];
                    if (ith_sample[feature_index] < node_value) {
                        // append to the left feature
                        split_res.m_feature_left.push_back(ith_sample);
                        // observation
                        split_res.m_obs_left.push_back(obs_y[loop_i]);
                    } else {
                        // append to the right
                        split_res.m_feature_right.push_back(ith_sample);
                        split_res.m_obs_right.push_back(obs_y[loop_i]);
                    }
                }
                
                // update terminal values
                if (m_type == AVERAGE) {
                    double mean_value = 0.0;
                    for (double obsL : split_res.m_obs_left) {
                        mean_value += obsL;
                    }
                    mean_value = mean_value / split_res.m_obs_left.size();
                    split_res.m_left_value = mean_value;
                    
                    mean_value = 0.0;
                    for (double obsR : split_res.m_obs_right) {
                        mean_value += obsR;
                    }
                    mean_value = mean_value / split_res.m_obs_right.size();
                    split_res.m_right_value = mean_value;
                    
                } else if (m_type == MAXIMAL) {
                    double max_value = 0.0;
                    VD::iterator iter = split_res.m_obs_left.begin();
                    if (++iter != split_res.m_obs_left.end()) {
                        max_value = *iter;
                    }
                    
                    while (++iter != split_res.m_obs_left.end()) {
                        double sel_value = *iter;
                        if (max_value < sel_value) {
                            max_value = sel_value;
                        }
                    }
                    
                    split_res.m_left_value = max_value;
                    
                    
                    // right value
                    max_value = 0.0;
                    iter = split_res.m_obs_right.begin();
                    if (++iter != split_res.m_obs_right.end()) {
                        max_value = *iter;
                    }
                    
                    while (++iter != split_res.m_obs_right.end()) {
                        double sel_value = *iter;
                        if (max_value < sel_value) {
                            max_value = sel_value;
                        }
                    }
                    
                    split_res.m_right_value = max_value;
                    
                } else {
                    // Unknown terminal type
                    assert(false);
                }
                
                // return the result
                return split_res;
            }
            
            /*
             *  The following function builds a regression tree from data
             */
            Node* buildTree(const VC<VD> &samples_x, const VD &obs_y) {
                
                // obtain the optimal split point
                m_current_depth = m_current_depth + 1;
                
                BestSplit best_split = findOptimalSplit(samples_x, obs_y);
                
                if (!best_split.m_status) {
                    if (m_current_depth > 0)
                        m_current_depth = m_current_depth - 1;
                    
                    return NULL;
                }
                
                // update feature importance info
                features_importance[best_split.m_feature_index] += 1;
                
                // split the data
                SplitRes split_data = splitData(samples_x, obs_y, best_split);
                
                // append current value to tree
                Node *new_node = new Node(best_split.m_node_value, best_split.m_feature_index, split_data.m_left_value, split_data.m_right_value);
                
                if (!m_root) {
                    m_root = new_node;
                    m_current_depth = 0;
                    // append left and right side
                    m_root->m_left_child = buildTree(split_data.m_feature_left, split_data.m_obs_left); // left
                    m_root->m_right_child = buildTree(split_data.m_feature_right, split_data.m_obs_right); // right
                } else {
                    // append left and right side
                    new_node->m_left_child = buildTree(split_data.m_feature_left, split_data.m_obs_left); // left
                    new_node->m_right_child = buildTree(split_data.m_feature_right, split_data.m_obs_right); // right
                }
                if (m_current_depth > 0)
                    m_current_depth--;
                
                return new_node;
            }
        };
        
        class PredictionForest {
        public:
            // class members
            double m_init_value;
            // the tree forest
            VC<RegressionTree> m_trees;
            // the learning rate
            double m_combine_weight;
            
            // construction function
            PredictionForest(double learning_rate) : m_init_value(0.0), m_combine_weight(learning_rate) {}
            
            /**
             * The method to make prediction for estimate of function's value from provided features
             *
             * @param feature_x the features to use for prediction
             * @return the estimated function's value
             */
            double predict(const VD &feature_x) {
                double re_res = m_init_value;
                
                if (m_trees.size() == 0) {
                    return re_res;
                }
                
                for (int i = 0; i < m_trees.size(); i++) {
                    re_res += m_combine_weight * m_trees[i].predict(feature_x);
                }
                
                return re_res;
            }
            
            /**
             * Calculates importance of each feature in input samples
             */
            VI featureImportances() {
                VI importances;
                for (int i = 0; i < m_trees.size(); i++) {
                    concatenate(importances, m_trees[i].features_importance);
                }
                return importances;
            }
        };
        
        
        class GradientBoostingMachine {
            // class members
            double m_sampling_size_ratio = 0.5;
            double m_learning_rate = 0.01;
            int m_tree_number = 100;
            
            // tree related parameters
            int m_tree_min_nodes = 10;
            int m_tree_depth = 3;
            
        public:
            
            GradientBoostingMachine(double sample_size_ratio, double learning_rate,
                                    int tree_number, int tree_min_nodes, int tree_depth) :
            m_sampling_size_ratio(sample_size_ratio), m_learning_rate(learning_rate), m_tree_number(tree_number),
            m_tree_min_nodes(tree_min_nodes), m_tree_depth(tree_depth) {
                // Check the validity of numbers
                Assert(sample_size_ratio > 0 && learning_rate > 0 && tree_number > 0 && tree_min_nodes >= 3 && tree_depth > 0,
                       "Wrong parameters");
            }
            
            /**
             * Fits a regression function using the Gradient Boosting Tree method.
             * On success, return function; otherwise, return null.
             *
             * @param input_x the input features
             * @param input_y the ground truth values - one per features row
             */
            PredictionForest *train(const VC<VD> &input_x, const VD &input_y) {
                
                // initialize forest
                PredictionForest *res_fun = new PredictionForest(m_learning_rate);
                
                // get the samples number
                size_t samples_num = input_y.size();
                
                Assert(samples_num == input_x.size() && samples_num > 0,
                       "Error: The input_x size should not be zero and should match the size of input_y");
                
                // holds indices of training data samples used for trees training
                // this will be used later for OOB error calculation
                VI used_indices(samples_num, -1);
                
                // get an initial guess of the function
                double mean_y = 0.0;
                for (double d : input_y) {
                    mean_y += d;
                }
                mean_y = mean_y / samples_num;
                res_fun->m_init_value = mean_y;
                
                
                // prepare the iteration
                VD h_value(samples_num);
                // initialize h_value
                int index = 0;
                while (index < samples_num) {
                    h_value[index] = mean_y;
                    index += 1;
                }
                
                // begin the boosting process
                int iter_index = 0;
                while (iter_index < m_tree_number) {
                    
                    // calculate the gradient
                    VD gradient;
                    index = 0;
                    for (double d : input_y) {
                        gradient.push_back(d - h_value[index]);
                        
                        // next
                        index++;
                    }
                    
                    // begin to sample
                    if (m_sampling_size_ratio < 0.99) {
                        // sample without replacement
                        
                        // we need to sample
                        RandomSample sampler((int)samples_num, (int) (m_sampling_size_ratio * samples_num));
                        
                        // get random index
                        VI sampled_index = sampler.get_sample_index();
                        
                        // data for growing trees
                        VC<VD> train_x;
                        VD train_y;
                        
                        for (int sel_index : sampled_index) {
                            // assign value
                            train_y.push_back(gradient[sel_index]);
                            train_x.push_back(input_x[sel_index]);
                            
                            // mark index as used
                            used_indices[sel_index] = 1;
                        }
                        
                        // fit a regression tree
                        RegressionTree tree;
                        
                        if (m_tree_depth > 0) {
                            tree.setDepth(m_tree_depth);
                        }
                        
                        if (m_tree_min_nodes > 0) {
                            tree.setMinNodes(m_tree_min_nodes);
                        }
                        
                        tree.buildRegressionTree(train_x, train_y);
                        
                        // store tree information
                        if (tree.m_root == NULL) {
                            // clear buffer
                            train_x.clear();
                            train_y.clear();
                            continue;
                        }
                        
                        res_fun->m_trees.push_back(tree);
                        
                        // update h_value information, prepare for the next iteration
                        int sel_index = 0;
                        while (sel_index < samples_num) {
                            h_value[sel_index] += m_learning_rate * tree.predict(input_x[sel_index]);
                            sel_index++;
                        }
                        
                    } else {
                        // use all data
                        // fit a regression tree
                        RegressionTree tree;
                        
                        // set parameters if needed
                        if (m_tree_depth > 0) {
                            tree.setDepth(m_tree_depth);
                        }
                        
                        if (m_tree_min_nodes > 0) {
                            tree.setMinNodes(m_tree_min_nodes);
                        }
                        
                        tree.buildRegressionTree(input_x, gradient);
                        
                        if (tree.m_root == NULL) {
                            // cannot update any more
                            break;
                        }
                        // store tree information
                        res_fun->m_trees.push_back(tree);
                        
                        // update h_value information, prepare for the next iteration
                        for (int loop_index = 0; loop_index < samples_num; loop_index++) {
                            h_value[loop_index] += m_learning_rate * tree.predict(input_x[loop_index]);
                        }
                    }
                    
                    // next iteration
                    iter_index++;
                }
                
                return res_fun;
            }
            
            PredictionForest *learnGradientBoostingRanker(const VC<VD> &input_x, const VC<VD> &input_y, const double tau) {
                PredictionForest *res_fun = new PredictionForest(m_learning_rate);
                
                size_t feature_num = input_x.size();
                
                Assert(feature_num == input_y.size() && feature_num > 0,
                       "The size of input_x should be the same as the size of input_y");
                
                VD h_value_x(feature_num, 0);
                VD h_value_y(feature_num, 0);
                
                int iter_index = 0;
                while (iter_index < m_tree_number) {
                    
                    // in the boosting ranker, randomly select half samples without replacement in each iteration
                    RandomSample sampler((int)feature_num, (int) (0.5 * feature_num));
                    
                    // get random index
                    VI sampled_index = sampler.get_sample_index();
                    
                    VC<VD> gradient_x;
                    VD gradient_y;
                    
                    for (int i = 0; i < sampled_index.size(); i++) {
                        int sel_index = sampled_index[i];
                        
                        gradient_x.push_back(input_x[sel_index]);
                        gradient_x.push_back(input_y[sel_index]);
                        
                        // get sample data
                        if (h_value_x[sel_index] < h_value_y[sel_index] + tau) {
                            double neg_gradient = h_value_y[sel_index] + tau - h_value_x[sel_index];
                            gradient_y.push_back(neg_gradient);
                            gradient_y.push_back(-1 * neg_gradient);
                        } else {
                            gradient_y.push_back(0.0);
                            gradient_y.push_back(0.0);
                        }
                        //                cerr << "sel_index: " << sel_index << endl;
                    }
                    
                    // fit a regression tree
                    RegressionTree tree;
                    //            tree.m_type = MAXIMAL;
                    
                    tree.buildRegressionTree(gradient_x, gradient_y);
                    
                    // store tree information
                    if (tree.m_root == NULL) {
                        continue;
                    }
                    
                    // update information
                    res_fun->m_trees.push_back(tree);
                    
                    double err = 0.0;
                    
                    for (int loop_index = 0; loop_index < feature_num; loop_index++) {
                        h_value_x[loop_index] += m_learning_rate * tree.predict(input_x[loop_index]);
                        h_value_y[loop_index] += m_learning_rate * tree.predict(input_y[loop_index]);
                        
                        if (h_value_x[loop_index] < h_value_y[loop_index] + tau) {
                            err += (h_value_x[loop_index] - h_value_y[loop_index] - tau) *
                            (h_value_x[loop_index] - h_value_y[loop_index] - tau);
                        }
                    }
                    
                    iter_index += 1;
                }
                
                
                
                return res_fun;
            }
        };
    }
}
#endif
