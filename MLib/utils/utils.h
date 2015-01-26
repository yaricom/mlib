//
//  Utils.h
//  activemoleculesC++
//
//  Created by Iaroslav Omelianenko on 1/9/15.
//  Copyright (c) 2015 yaric. All rights reserved.
//

#ifndef activemoleculesC___Utils_h
#define activemoleculesC___Utils_h

#include "defines.h"
#include "matrix.h"

using namespace std;

static bool LOG_DEBUG = true;
/*! the message buffer length */
const int kPrintBuffer = 1 << 12;

inline std::string spf(const char *fmt, ...) {
    std::string msg(kPrintBuffer, '\0');
    va_list args;
    va_start(args, fmt);
    vsnprintf(&msg[0], kPrintBuffer, fmt, args);
    va_end(args);
    
    return msg;
}

inline void Printf(const char *fmt, ...) {
    if (LOG_DEBUG) {
        std::string msg(kPrintBuffer, '\0');
        va_list args;
        va_start(args, fmt);
        vsnprintf(&msg[0], kPrintBuffer, fmt, args);
        va_end(args);
        fprintf(stderr, "%s", msg.c_str());
    }
}

inline void Assert(bool exp, const char *fmt, ...) {
    if (!exp) {
        std::string msg(kPrintBuffer, '\0');
        va_list args;
        va_start(args, fmt);
        vsnprintf(&msg[0], kPrintBuffer, fmt, args);
        va_end(args);
        fprintf(stderr, "AssertError:%s\n", msg.c_str());
        exit(-1);
    }
}



// toString
template<class T> std::string i2s(T x) {std::ostringstream o; o << x; return o.str();}
// print vector
template<class T> void print(std::vector < T > v) {cerr << "[";if (v.size()) cerr << v[0];FOR(i, 1, v.size()) cerr << ", " << v[i];cerr << "]" << endl;}

// print matrix
inline void print(Matrix &m) {
    size_t rows = m.rows();
    size_t cols = m.cols();
    cerr << "[" << endl;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cerr << m(i, j);
            if (j < cols - 1) {
                cerr << ",\t";
            }
        }
        cerr << endl;
    }
    cerr << "]" << endl;
}

// print vector
inline void print(Vector &m) {
    size_t size = m.size();
    cerr << "[" << endl;
    for (int i = 0; i < size; i++) {
        cerr << m[i];
        if (i < size - 1) {
            cerr << ",\t";
        }
        
        cerr << endl;
    }
    cerr << "]" << endl;
}


#endif
