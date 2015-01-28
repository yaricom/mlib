//
//  io.h
//  MLib
//
//  Created by Iaroslav Omelianenko on 1/28/15.
//  Copyright (c) 2015 yaric. All rights reserved.
//

#ifndef MLib_io_h
#define MLib_io_h

#include <stdio.h>
#include <exception>
#include "matrix.h"

namespace nologin {
    namespace utils {
        
        using namespace math;
        /**
         * Writes specified matrix in LibSVM format.
         * @param fileName The path to the file to be written
         * @param mat The matrix to be written
         * @param classCol The column in matrix holding class values.(Default: last one)
         */
        void storeMatrixAsLibSVM(const char* fileName, const Matrix &mat, int classCol = -1) {
            FILE *fp;
            if (!(fp = fopen(fileName, "w"))) {
                throw runtime_error("Failed to open file!");
            }
            
            if (classCol < 0) {
                classCol = (int)mat.cols() - 1;
            }
            assert(classCol < mat.cols());
            // write to the buffer
            for (int row = 0; row < mat.rows(); row++) {
                // write class value first
                double val = mat(row, classCol);
                fprintf(fp, "%f", val);
                for (int col = 1; col < mat.cols(); col++) {
                    val = mat(row, col);
                    if (val) {
                        // write only non zero
                        fprintf(fp, " %d:%f", col, val);
                    }
                }
                fprintf(fp, "\n");
            }
            
            // close file
            fclose(fp);
        }
        
        char *strsep(char **stringp, const char *delim){
            char *s;
            const char *spanp;
            int c, sc;
            char *tok;
            
            if ((s = *stringp) == NULL)
                return (NULL);
            for (tok = s;;) {
                c = *s++;
                spanp = delim;
                do {
                    if ((sc = *spanp++) == c) {
                        if (c == 0)
                            s = NULL;
                        else
                            s[-1] = 0;
                        *stringp = s;
                        return (tok);
                    }
                } while (sc != 0);
            }
        }
        
        /**
         * Loads matrix from libsvm file
         * @param fileName The path to the file with matrix data
         * @param classCol The column in matrix holding class values.(Default: last one)
         */
        Matrix& loadMatrixFromLibSVM(const char* fileName, int classCol = -1) {
            FILE *fp;
            if (!(fp = fopen(fileName, "r"))) {
                throw runtime_error("Failed to open file!");
            }
            VDD data;
            VD classVal;
            char buf[2048];
            while (fscanf(fp,"%[^\n]\n", buf) == 1) {
                int i = 0, index = 0;
                float val = 0;
                char *ptr=buf;
                VD row;
                while (ptr != NULL) {
                    char *token = strsep(&ptr, " ");
                    if (strlen(token) != 0) {
                        if (index++ == 0) {
                            // read class value
                            classVal.push_back(atof(token));
                        } else {
                            // read feature value
                            sscanf(token, "%d:%f", &i, &val);
                            if (i > row.size()) {
                                row.resize(i, 0);
                            }
                            row[i - 1] = val;
                        }
                    }
                }
                
                // store row
                data.push_back(row);
            }
            
#warning            // pass through the collected data and ckeck dimensions and append class values
            
            
            // close file
            fclose(fp);
            
            return *(new Matrix(data));
        }
    }
}

#endif
