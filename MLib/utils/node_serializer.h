//
//  node_serializer.h
//  MLib
//
//  Created by Iaroslav Omelianenko on 8/6/15.
//  Copyright (c) 2015 yaric. All rights reserved.
//

#ifndef MLib_node_serializer_h
#define MLib_node_serializer_h

#include <ostream>
#include <iostream>
#include <vector>

using namespace std;

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

void serialize(const Node *node, ostream &os) {
    if (node == NULL) {
        os << ".";
        return;
    }
    // store node values
    os << node->m_node_value;
    os << " ";
    os << node->m_feature_index;
    os << " ";
    os << node->m_terminal_left;
    os << " ";
    os << node->m_terminal_right;
    os << " ";
    
    // proceed with childs
    serialize(node->m_left_child, os);
    serialize(node->m_right_child, os);
}

Node *deserialize(istream &is) {
    char c = is.get();
    if (c == '.') {
        return NULL;
    }
    // return charcter back
    is.unget();
    // read next node
    double node_value;
    is >> node_value;
    is.get();
    int feature_index;
    is >> feature_index;
    is.get();
    double terminal_left;
    is >> terminal_left;
    is.get();
    double terminal_right;
    is >> terminal_right;
    is.get();
    Node *node = new Node(node_value, feature_index, terminal_left, terminal_right);
    
    // proceed with childs
    node->m_left_child = deserialize(is);
    node->m_right_child = deserialize(is);
    
    return node;
}

#endif
