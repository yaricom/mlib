//
//  Matrix.h
//  activemoleculesC++
//
//  Created by Iaroslav Omelianenko on 1/20/15.
//  Copyright (c) 2015 yaric. All rights reserved.
//

#ifndef activemoleculesC___Matrix_h
#define activemoleculesC___Matrix_h

#include "defines.h"

using namespace std;

/*
 * The simple matrix implementation
 */
class Matrix {
    
protected:
    /** Array for internal storage of elements. */
    VDD A;
    
    /** Row and column dimensions.
     * m row dimension.
     * n column dimension.
     */
    size_t m, n;
    
    
public:
    
    /**
     * Construct an m-by-n matrix of zeros.
     *
     * @param rows Number of rows.
     * @param cols Number of colums.
     */
    
    Matrix(const size_t rows, const size_t cols) {
        m = rows;
        n = cols;
        for (int i = 0; i < m; i++) {
            VD row(n, 0);
            A.push_back(row);
        }
    }
    
    /**
     * Construct a matrix from a 2-D array.
     *
     * @param A Two-dimensional array of doubles.
     */
    Matrix(const VDD &arr) {
        m = arr.size();
        n = arr[0].size();
        for (int i = 0; i < m; i++) {
            assert(arr[i].size() == n);
        }
        A = arr;
    }
    
    /**
     * Construct a matrix from a one-dimensional packed array
     *
     * @param vals One-dimensional array of doubles, packed by columns (ala Fortran).
     * @param rows Number of rows.
     */
    Matrix(const VD &vals, const size_t rows) {
        m = rows;
        n = (m != 0 ? vals.size() / m : 0);
        assert (m * n == vals.size());
        
        for (int i = 0; i < m; i++) {
            VD row(n);
            for (int j = 0; j < n; j++) {
                row[j] = vals[i + j * m];
            }
            A.push_back(row);
        }
    }
    
    /**
     * Construct a matrix from a one-dimensional packed array
     *
     * @param vals One-dimensional array of doubles, packed by rows.
     * @param rows Number of rows.
     */
    Matrix(const size_t cols, const VD &vals) {
        n = cols;
        m = (n != 0 ? vals.size() / n : 0);
        assert (m * n == vals.size());
        
        for (int i = 0; i < m; i++) {
            VD row(n);
            for (int j = 0; j < n; j++) {
                row[j] = vals[i * n + j];
            }
            A.push_back(row);
        }
        
    }
    
    /**
     * Get row dimension.
     *
     * @return m, the number of rows.
     */
    size_t rows() const {
        return m;
    }
    
    /**
     * Get column dimension.
     *
     * @return n, the number of columns.
     */
    size_t cols() const {
        return n;
    }
    
    /**
     * Make a one-dimensional column packed copy of the internal array.
     *
     * @return Matrix elements packed in a one-dimensional array by columns.
     */
    void columnPackedCopy(VD &vals) {
        vals.resize(m * n, 0);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                vals[i + j * m] = A[i][j];
            }
        }
    }
    
    /**
     * Make a one-dimensional row packed copy of the internal array.
     *
     * @return Matrix elements packed in a one-dimensional array by rows.
     */
    void rowPackedCopy(VD &vals) {
        vals.resize(m * n, 0);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                vals[i * n + j] = A[i][j];
            }
        }
    }
    
    /**
     * Adds specified row to the end of matrix
     */
    void addRow(const VD &row) {
        assert(row.size() == n);
        A.push_back(row);
    }
    
    /**
     * Get a submatrix.
     *
     * @param i0 Initial row index
     * @param i1 Final row index
     * @param j0 Initial column index
     * @param j1 Final column index
     * @return A(i0:i1,j0:j1)
     */
    
    Matrix subMatrix(int i0, int i1, int j0, int j1) {
        assert(i0 >= 0 && i0 < i1 && i1 < m && j0 >= 0 && j0 < j1 && j1 < n);
        Matrix X(i1 - i0 + 1, j1 - j0 + 1);
        for (int i = i0; i <= i1; i++) {
            for (int j = j0; j <= j1; j++) {
                X(i - i0, j - j0) = A[i][j];
            }
        }
        return X;
    }
    
    /**
     * Get a submatrix.
     *
     * @param i0 Initial row index
     * @param i1 Final row index
     * @param c Array of column indices.
     * @return A(i0:i1,c(:))
     */
    Matrix subMatrix(const int i0, const int i1, const VI &c) {
        assert(i0 >= 0 && i0 < i1 && i1 < m);
        Matrix X(i1 - i0 + 1, c.size());
        for (int i = i0; i <= i1; i++) {
            for (int j = 0; j < c.size(); j++) {
                assert(c[j] < n && c[j] >= 0);
                X(i - i0, j) = A[i][c[j]];
            }
        }
        return X;
    }
    
    /**
     * Get a submatrix.
     *
     * @param r
     *            Array of row indices.
     * @param j0
     *            Initial column index
     * @param j1
     *            Final column index
     * @return A(r(:),j0:j1)
     * @exception ArrayIndexOutOfBoundsException
     *                Submatrix indices
     */
    
    Matrix subMatrix(const VI &r, const int j0, const int j1) {
        assert(j0 >= 0 && j0 < j1 && j1 < n);
        Matrix X(r.size(), j1 - j0 + 1);
        for (int i = 0; i < r.size(); i++) {
            assert(r[i] < m && r[i] >= 0);
            for (int j = j0; j <= j1; j++) {
                X(i, j - j0) = A[r[i]][j];
            }
        }
        return X;
    }
    
    /**
     * Set a submatrix.
     *
     * @param i0 Initial row index
     * @param i1 Final row index
     * @param j0 Initial column index
     * @param j1 Final column index
     * @param X A(i0:i1,j0:j1)
     */
    void setMatrix(const int i0, const int i1, const int j0, const int j1, const Matrix &X) {
        assert(i0 >= 0 && i0 < i1 && i1 < m && j0 >= 0 && j0 < j1 && j1 < n);
        for (int i = i0; i <= i1; i++) {
            for (int j = j0; j <= j1; j++) {
                A[i][j] = X.A[i - i0][j - j0];
            }
        }
    }
    
    /**
     * Set a submatrix.
     *
     * @param r Array of row indices.
     * @param c Array of column indices.
     * @param X A(r(:),c(:))
     */
    void setMatrix(const VI &r, const VI &c, const Matrix &X) {
        for (int i = 0; i < r.size(); i++) {
            assert(r[i] < m && r[i] >= 0);
            for (int j = 0; j < c.size(); j++) {
                assert(c[j] < n && c[j] >= 0);
                A[r[i]][c[j]] = X.A[i][j];
            }
        }
    }
    
    /**
     * Set a submatrix.
     *
     * @param r Array of row indices.
     * @param j0 Initial column index
     * @param j1 Final column index
     * @param X A(r(:),j0:j1)
     */
    void setMatrix(const VI &r, const int j0, const int j1, const Matrix &X) {
        assert(j0 >=0 && j0 < j1 && j1 < n);
        for (int i = 0; i < r.size(); i++) {
            assert(r[i] < m && r[i] >= 0);
            for (int j = j0; j <= j1; j++) {
                A[r[i]][j] = X.A[i][j - j0];
            }
        }
    }
    
    /**
     * Set a submatrix.
     *
     * @param i0
     *            Initial row index
     * @param i1
     *            Final row index
     * @param c
     *            Array of column indices.
     * @param X
     *            A(i0:i1,c(:))
     * @exception ArrayIndexOutOfBoundsException
     *                Submatrix indices
     */
    
    void setMatrix(const int i0, const int i1, const VI c, const Matrix &X) {
        assert(i0 >= 0 && i0 < i1 && i1 < m);
        for (int i = i0; i <= i1; i++) {
            for (int j = 0; j < c.size(); j++) {
                assert(c[j] < n && c[j] >= 0);
                A[i][c[j]] = X.A[i - i0][j];
            }
        }
    }
    
    /**
     * Matrix transpose.
     *
     * @return A'
     */
    Matrix transpose() {
        Matrix X(n, m);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                X(j, i) = A[i][j];
            }
        }
        return X;
    }
    
#pragma mark - Algebraic functions
    /**
     * Scales this matrix to have all samles scaled to fit range [min, max]
     *
     * X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
     * X_scaled = X_std * (max - min) + min
     */
    void scaleMinMax(double min, double max) {
        // fin min/max per sample per feature
        VD fMins(n, numeric_limits<double>().max()), fMaxs(n, numeric_limits<double>().min());
        for (int i = 0 ; i < m; i++) {
            for (int j = 0; j < n; j++) {
                fMaxs[j] = std::max<double>(fMaxs[j], A[i][j]);
                fMins[j] = std::min<double>(fMins[j], A[i][j]);
            }
        }
        
        // find X scaled
        for (int i = 0 ; i < m; i++) {
            for (int j = 0; j < n; j++) {
                double X = A[i][j], X_min = fMins[j], X_max = fMaxs[j];
                double X_std = (X - X_min) / (X_max - X_min);
                A[i][j] = X_std * (max - min) + min;
            }
        }
    }
    
    /**
     * One norm
     *
     * @return maximum column sum.
     */
    double norm1() const {
        double f = 0;
        for (int j = 0; j < n; j++) {
            double s = 0;
            for (int i = 0; i < m; i++) {
                s += abs(A[i][j]);
            }
            f = max(f, s);
        }
        return f;
    }
    
    /**
     * Infinity norm
     *
     * @return maximum row sum.
     */
    double normInf() const {
        double f = 0;
        for (int i = 0; i < m; i++) {
            double s = 0;
            for (int j = 0; j < n; j++) {
                s += abs(A[i][j]);
            }
            f = max(f, s);
        }
        return f;
    }
    
    /**
     * Frobenius norm
     *
     * @return sqrt of sum of squares of all elements.
     */
    double normF() const {
        double f = 0;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                f = hypot(f, A[i][j]);
            }
        }
        return f;
    }
    
    /**
     * Unary minus
     *
     * @return -A
     */
    Matrix uminus() const {
        Matrix X(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                X(i, j) = -A[i][j];
            }
        }
        return X;
    }
 
#pragma mark - Operators
    /**
     * Copy matrix B which should have equal dimensions. A = B
     */
    Matrix& operator=(const Matrix &B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                this->A[i][j] = B.A[i][j];
            }
        }
        return *this;
    }
    
    /**
     * Returns reference to element at A(i,j)
     * @param i Row index.
     * @param j Column index.
     */
    double& operator()(const int i, const int j) {
        assert( i >= 0 && i < m && j >=0 && j < n);
        return A[i][j];
    }
    
    /**
     * Returns reference to element at A(i,j)
     * @param i Row index.
     * @param j Column index.
     */
    double operator()(const int i, const int j) const{
        assert( i >= 0 && i < m && j >=0 && j < n);
        return A[i][j];
    }
    
    /**
     * C = A + B
     *
     * @param B another matrix
     * @return A + B
     */
    Matrix operator+(const Matrix& B) const {
        checkMatrixDimensions(B);
        Matrix X(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                X(i, j) = this->A[i][j] + B.A[i][j];
            }
        }
        return X;
    }
    
    /**
     * A = A + B
     *
     * @param B
     *            another matrix
     * @return A + B
     */
    
    Matrix& operator+=(const Matrix &B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                this->A[i][j] = this->A[i][j] + B.A[i][j];
            }
        }
        return *this;
    }
    
    /**
     * C = A - B
     *
     * @param B another matrix
     * @return A - B
     */
    
    Matrix operator-(const Matrix &B) const {
        checkMatrixDimensions(B);
        Matrix X(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                X(i, j) = this->A[i][j] - B.A[i][j];
            }
        }
        return X;
    }
    
    /**
     * A = A - B
     *
     * @param B
     *            another matrix
     * @return A - B
     */
    
    Matrix& operator-=(const Matrix &B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                this->A[i][j] = this->A[i][j] - B.A[i][j];
            }
        }
        return *this;
    }
    
    /**
     * Element-by-element multiplication, C = A.*B
     *
     * @param B another matrix
     * @return A.*B
     */
    Matrix operator*(const Matrix &B) const {
        checkMatrixDimensions(B);
        Matrix X(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                X(i, j) = this->A[i][j] * B.A[i][j];
            }
        }
        return X;
    }
    
    /**
     * Element-by-element multiplication in place, A = A.*B
     *
     * @param B another matrix
     * @return A.*B
     */
    
    Matrix& operator*=(const Matrix &B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                this->A[i][j] = this->A[i][j] * B.A[i][j];
            }
        }
        return *this;
    }
    
    /**
     * Multiply a matrix by a scalar, C = s*A
     *
     * @param s scalar
     * @return s*A
     */
    
    Matrix operator*(const double s) const {
        Matrix X(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                X(i, j) = s * this->A[i][j];
            }
        }
        return X;
    }
    
    /**
     * Multiply a matrix by a scalar in place, A = s*A
     *
     * @param s scalar
     * @return replace A by s*A
     */
    
    Matrix& operator*=(const double s) {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                this->A[i][j] = s * this->A[i][j];
            }
        }
        return *this;
    }
    
    /**
     * Element-by-element right division, C = A./B
     *
     * @param B
     *            another matrix
     * @return A./B
     */
    
    Matrix operator/(const Matrix &B) const {
        checkMatrixDimensions(B);
        Matrix X(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                X(i, j) = this->A[i][j] / B.A[i][j];
            }
        }
        return X;
    }
    
    /**
     * Element-by-element right division in place, A = A./B
     *
     * @param B
     *            another matrix
     * @return A./B
     */
    
    Matrix& operator/=(const Matrix &B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                this->A[i][j] = this->A[i][j] / B.A[i][j];
            }
        }
        return *this;
    }
    
    /**
     * Element-by-element right division by scalar in place, A = A./s
     *
     * @param B
     *            another matrix
     * @return A./B
     */
    
    Matrix& operator/=(const double s) {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                this->A[i][j] =  this->A[i][j] / s;
            }
        }
        return *this;
    }
    
    /**
     * Division of matrix by a scalar, C = A/s
     *
     * @param s scalar
     * @return s*A
     */
    Matrix operator/(const double s) const {
        Matrix X(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                X(i, j) = this->A[i][j] / s;
            }
        }
        return X;
    }
    
    /**
     * Element-by-element equality check
     */
    bool operator==(const Matrix &B) const {
        if (m != B.m && n != B.n) {
            return false;
        }
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                if (A[i][j] != B.A[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }
    
    /**
     * Element-by-element left division, C = A.\B
     *
     * @param B another matrix
     * @return A.\B
     */
    Matrix arrayLeftDivide(const Matrix &B) const {
        checkMatrixDimensions(B);
        Matrix X(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                X(i, j) = B.A[i][j] / this->A[i][j];
            }
        }
        return X;
    }
    
    /**
     * Element-by-element left division in place, A = A.\B
     *
     * @param B
     *            another matrix
     * @return A.\B
     */
    
    Matrix& arrayLeftDivideEquals(const Matrix &B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                this->A[i][j] = B.A[i][j] / this->A[i][j];
            }
        }
        return *this;
    }
    
    /**
     * Linear algebraic matrix multiplication, A * B
     * Matrix product
     *
     * @param B another matrix
     * @return Matrix product, A x B
     */
    
    Matrix matmul(const Matrix &B) const {
        // Matrix inner dimensions must agree.
        assert (B.m == n);
        
        Matrix X(m, B.n);
        double Bcolj[n];
        for (int j = 0; j < B.n; j++) {
            for (int k = 0; k < n; k++) {
                Bcolj[k] = B.A[k][j];
            }
            for (int i = 0; i < m; i++) {
                VD Arowi = A[i];
                double s = 0;
                for (int k = 0; k < n; k++) {
                    s += Arowi[k] * Bcolj[k];
                }
                X(i, j) = s;
            }
        }
        return X;
    }
    
    
    
protected:
    void checkMatrixDimensions(Matrix B) const {
        assert (B.m != m || B.n != n);
    }
};

/**
 * The simple Vector implementation
 */
class Vector : Matrix {
public:
    /**
     * Constructs zero vector with the specified size.
     * @param size the capacity of the vector
     */
    Vector(int size) : Matrix(size, 1) {}
    
    /**
     * Constructs a vector with the specified list of values.
     * @param lst a list of values
     */
    Vector(VD lst) : Matrix(lst.size(), 1) {
        for (int i = 0; i < lst.size(); i++) {
            A[i][0] = lst[i];
        }
    }
    
    /**
     * Returns the number of values in this vector.
     * @return the number of values in this vector.
     */
    size_t size() const {
        return rows();
    }
    
#pragma mark - Operators
    /**
     * Returns reference to element at A(i,j)
     * @param i the index.
     */
    double& operator()(const int i) {
        assert( i >= 0 && i < m);
        return A[i][0];
    }
    
    /**
     * Copy vector B which should have equal dimensions. A = B
     */
    Vector& operator=(const Vector &B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < m; i++) {
            this->A[i][0] = B.A[i][0];
        }
        return *this;
    }
    
    /**
     * Returns the result of vector subtraction.
     *
     * @param vthe vector to subtract
     * @return the result of vector subtraction.
     */
    Vector operator-(const Vector &v) {
        checkMatrixDimensions(v);
        int size = (int)this->size();
        Vector result(size);
        for (int i = 0; i < size; i++) {
            result(i) = result(i) - v.A[i][0];
        }
        return result;
    }
    
    /**
     * Returns the result of vector addition.
     *
     * @param v
     *            the vector to add
     *
     * @return the result of vector addition.
     */
    Vector operator+(const Vector &v) {
        checkMatrixDimensions(v);
        int size = (int)this->size();
        Vector result(size);
        for (int i = 0; i < size; i++) {
            result(i) = result(i) + v.A[i][0];
        }
        return result;
    }
    
    /**
     * Element-by-element right division by scalar in place, A = A./s
     *
     * @param B another matrix
     * @return A./B
     */
    
    Vector& operator/=(const double s) {
        for (int i = 0; i < m; i++) {
            this->A[i][0] = this->A[i][0] / s;
        }
        return *this;
    }
    /**
     * Element-by-element equality check
     */
    bool operator==(const Vector &B) const {
        if (this->size() != B.size()) {
            return false;
        }
        for (int i = 0; i < m; i++) {
            if (this->A[i][0] != B.A[i][0]) {
                return false;
            }
        }
        return true;
    }
    
};

#pragma mark - Matrix functions
/**
 * Calculates matrix mean by columns
 */
Vector& mean(const Matrix &m){
    size_t cols = m.cols(), rows = m.rows();
    Vector *current = new Vector((int)cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            (*current)(j) += m(i, j);
        }
    }
    (*current) /= rows;
    return *current;
}

/**
 * Calculates variance of matrix by columns.
 *
 * @return the row vector containing the variance of the elements per column.
 */
Vector& variance(const Matrix &m){
    Vector vmean = mean(m);
    size_t cols = m.cols(), rows = m.rows();
    
    Vector *X = new Vector((int)cols);
    double diff;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            diff = m(i, j) - vmean(j);
            (*X)(j) += diff * diff;
        }
    }

    return *X;
}

/**
 * Calculates standard deviation of matrix by columns.
 *
 * @return the row vector containing the standard deviation of the elements of each column.
 */
Vector& stdev(const Matrix &m)  {
    Vector &var = variance(m);
    size_t n = var.size();

    for (int j = 0; j < n; j++) {
        var(j) = sqrt(var(j) / (n - 1));
    }
    return var;
}

#endif
