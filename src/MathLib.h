#ifndef __MATHLIB_H__
#define __MATHLIB_H__

#pragma once

#include <cmath>
#include <math.h>
#include <random>
#include <cassert>
#include <iostream>
#include <algorithm>

static const double INF = std::numeric_limits<double>::infinity();
static const double EPS = std::numeric_limits<double>::epsilon();
static const double PI  = 3.141592653589793115997963468544185161590576171875;

/********************************************************************************
*                               Vector Defination                               *
********************************************************************************/
// template parameter used here
template<int n, typename T> struct vec {
    T data[n] = {0};

    vec() = default;
    T & operator[] (const int i)       { assert(i >= 0 && i < n); return data[i]; }
    T   operator[] (const int i) const { assert(i >= 0 && i < n); return data[i]; }
    inline vec operator-() const { 
        for (int i = 0; i < n; data[i] = -data[i])
        return *this;
    }
    inline vec& operator+=(const vec<n,T> &v) {
        for (int i = 0; i < n; data[i] += v[i])
        return *this;
    }
    inline vec& operator*=(T t) {
        for (int i = 0; i < n; data[i] *= t)
        return *this;
    }
    inline vec& operator/=(T t) { return *this *= 1/t; }

    inline T norm2() const { return *this * *this; }
    inline T norm()  const { return std::sqrt(norm2()); }
    inline vec normalize() const { return (*this) * (1 / norm()); }
};

// overload operator <<, vector print
template<int n,typename T> std::ostream& operator<< (std::ostream& out, const vec<n,T>& v) {
    std::cout << "[ ";
    for (int i = 0; i < n; i++) std::cout << v[i] << " ";
    std::cout << ']';
    return out;
};

// overload operator +, vector - vector addition
template<int n, typename T>
inline vec<n,T> operator+ (const vec<n,T>& u, const vec<n,T>& v) {
    vec<n,T> res = u;
    for (int i = n; i--; res[i] += v[i]) ;
    return res;
};

// overload operator +, vector - scalar addition
template<int n, typename T>
inline vec<n,T> operator+ (const vec<n,T>& u, const T& k) {
    vec<n,T> res = u;
    for (int i = n; i--; res[i] += k) ;
    return res;
};

// overload operator -, vector - vector subtraction
template<int n, typename T>
inline vec<n,T> operator- (const vec<n,T>& u, const vec<n,T>& v) {
    vec<n,T> res = u;
    for (int i = n; i--; res[i] -= v[i]) ;
    return res;
};

// overload operator -, vector - scalar subtraction
template<int n, typename T>
inline vec<n,T> operator- (const vec<n,T>& u, const T x) {
    vec<n,T> res = u;
    for (int i = n; i--; res[i] -= x) ;
    return res;
};

// overload operator *, vector dot product
template<int n, typename T>
inline T operator* (const vec<n,T>& u, const vec<n,T>& v) {
    T res = 0;
    for (int i = n; i--; res += u[i] * v[i]) ;
    return res;
};

// overload operator *, vector - scalar multipication
template<int n, typename T>
inline vec<n,T> operator* (const T& k, const vec<n,T>& u) {
    vec<n,T> res = u;
    for (int i = n; i--; res[i] *= k) ;
    return res;
};

// overload operator *, vector - scalar multipication
template<int n, typename T>
inline vec<n,T> operator* (const vec<n,T>& u, const T& k) {
    return k * u;
};

// overload operator /, vector - scalar division
template<int n, typename T> vec<n,T>
inline operator/ (const vec<n,T>& u, const T& k) {
    return (1 / k) * u;
};

// vector element-wise product
template<int n, typename T> vec<n,T>
inline mult(const vec<n,T> &v1, const vec<n,T> &v2) {
    vec<n,T> res;
    for (int i = n; i--; res[i] = v1[i] * v2[i]) ;
    return res;
}

template<int n1, int n2, typename T>
inline vec<n1,T> embed(const vec<n2,T>& v, T fill = 0) {
    vec<n1,T> res;
    for (int i = n1; i--; res[i] = (i < n2) ? v[i] : fill) ;
    return res;
}

template<int n1, int n2, typename T>
inline vec<n1,T> slice(const vec<n2,T> v, const int m=0, const int n=n1-1) {
    assert(n2 > (n - m + 1) && n >= m);
    vec<n1,T> res;
    for (int i = 0, j = m; j <= n; i++, j++) res[i] = v[j];
    return res;
}

// compute the projection of u onto v
template<int n, typename T> 
inline vec<n,T> proj(const vec<n,T>& u, const vec<n,T>& v) {
    vec<n,T> res = u;
    vec<n,T> vNorm = (1 / v.norm()) * v;
    return (res *  vNorm) * vNorm;
}

// Specialization for vec<2,double>
template<> struct vec<2,double> {
    double x{}, y{};

    vec() = default;
    vec(double x_, double y_): x(x_), y(y_) {}
    double & operator[] (const int i)       { assert(i >= 0 && i < 2); return i ? y : x; }
    double   operator[] (const int i) const { assert(i >= 0 && i < 2); return i ? y : x; }
    inline vec operator-() const { return {-x, -y}; }
    inline vec& operator+=(const vec<2,double> &v) {
        x += v.x;
        y += v.y;
        return *this;
    }
    inline vec& operator*=(double t) {
        x *= t;
        y *= t;
        return *this;
    }
    inline vec& operator/=(double t) { return *this *= 1/t; }

    inline double norm2() const { return x*x + y*y; }
    inline double norm()  const { return sqrt(norm2()); }
    inline vec normalize() const { return (*this) * (1 / norm()); }
};

// Specialization for vec<3,double>
template<> struct vec<3,double> {
    double x{}, y{}, z{};

    vec() = default;
    vec(double x_, double y_, double z_ = 0.0): x(x_), y(y_), z(z_) {}
    double & operator[] (const int i)       { assert(i >= 0 && i < 3); return i ? (i == 2 ? z : y) : x; }
    double   operator[] (const int i) const { assert(i >= 0 && i < 3); return i ? (i == 2 ? z : y) : x; }
    inline vec  operator-() const { return {-x, -y, -z}; }
    inline vec& operator+=(const vec &v) {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    inline vec& operator*=(double t) {
        x *= t;
        y *= t;
        z *= t;
        return *this;
    }
    inline vec& operator/=(double t) { return *this *= 1/t; }
    inline bool operator<(const vec<3,double> &v) const {
        return ((x < v.x) && (y < v.y) && (z < v.z));
    }

    inline double norm2() const { return x*x + y*y + z*z; }
    inline double norm()  const { return sqrt(norm2()); }
    inline vec normalize() const { return (*this) * (1 / norm()); }
    // vector compnent-wise product
    inline vec mult(const vec<3,double> &v) {
        return {x*v.x, y*v.y, z*v.z};
    }
    // vector cross product
    inline vec cross(const vec<3,double> &v) const {
        return {y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x};
    }
};

// overload operator +, vec3 + vec3
inline vec<3,double> operator+(const vec<3,double> &u, const vec<3,double> &v) {
    return {u.x + v.x, u.y + v.y, u.z + v.z};
}

// overload operator +, vec3 + double
inline vec<3,double> operator+(const vec<3,double> &u, const double &k) {
    return {u.x + k, u.y + k, u.z + k};
}

// overload operator -, vec3 - vec3
inline vec<3,double> operator-(const vec<3,double> &u, const vec<3,double> &v) {
    return {u.x - v.x, u.y - v.y, u.z - v.z};
}

// overload operator -, vec3 - double
inline vec<3,double> operator-(const vec<3,double> &u, const double &k) {
    return {u.x - k, u.y - k, u.z - k};
}

// overload operator *, double * vec3
inline vec<3,double> operator*(const double &k, const vec<3,double> &u) {
    return {u.x * k, u.y * k, u.z * k};
}

// overload operator *, vec3 * double addition
inline vec<3,double> operator*(const vec<3,double> &u, const double &k) {
    return k * u;
}

// overload operator *, vec3 dot vec3
inline double operator*(const vec<3,double> &u, const vec<3,double> &v) {
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

// overload operator /, vec3 / double
inline vec<3,double> operator/(const vec<3,double> &v, const double &k) {
    return (1/k) * v;
}



// vector compnent-wise product
inline vec<3,double> mult(const vec<3,double> &v1, const vec<3,double> &v2) {
    return {v1.x*v2.x, v1.y*v2.y, v1.z*v2.z};
}

// vector cross product
inline vec<3,double> cross(const vec<3,double> &v1, const vec<3,double> &v2) {
    return {v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x};
}

// Specialization for vec<4,double>
template<> struct vec<4,double> {
    double x{}, y{}, z{}, w{};

    vec() = default;
    vec(double x_, double y_, double z_, double w_ = 0.0): x(x_), y(y_), z(z_), w(w_) {}
    vec(vec<3,double> v, double w=0.0): x(v.x), y(v.y), z(v.z), w(w) {}
    double & operator[] (const int i)       { assert(i >= 0 && i < 4); return i ? (i == 3 ? w : (i == 2 ? z : y)) : x; }
    double   operator[] (const int i) const { assert(i >= 0 && i < 4); return i ? (i == 3 ? w : (i == 2 ? z : y)) : x; }

    inline vec<3,double> xyz() { return {x, y, z}; }
    inline double norm2() const { return x*x + y*y + z*z + w*w; }
    inline double norm()  const { return sqrt(norm2()); }
    inline vec normalize() const { return (*this) * (1 / norm()); }
};

inline double operator*(const vec<4,double> &u, const vec<4,double> &v) {
    return u.x * v.x + u.y * v.y + u.z * v.z + u.w * v.w;
}

typedef vec<2,double> vec2;
typedef vec<3,double> vec3;
typedef vec<4,double> vec4;

/********************************************************************************
*                               Matrix Defination                               *
********************************************************************************/
template<int nrows, int ncols, typename T> struct mat;

template<typename T> T detRecursive(mat<2,2,T> m) {
    return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

template<int N, typename T> T detRecursive(mat<N,N,T> m) {
    mat<N-1,N-1,T> submatrix;
    vec<N-1,T> row;
    T det = 0.0;

    for (int j = 0; j < N; ++j) {
        for (int i = 1; i < N; ++i) {
            for (int k = 0; k < N; ++k) {
                if (k != j) {
                    row[i-1] = m[i][k];
                }
            }
            submatrix[i-1] = row;
        }

        det += m[0][j] * detRecursive(submatrix) * ((j % 2 == 0) ? 1 : -1);
    }

    return det;
}

template<int N, typename T> T detGauss(mat<N,N,T> m) {
    int swaps = 0;

    for (int i = 0; i < N - 1; ++i) {
        int pivot_row = i;
        for (int j = i + 1; j < N; ++j) {
            if (std::abs(m[j][i]) > std::abs(m[pivot_row][i])) {
                pivot_row = j;
            }
        }

        if (pivot_row != i) {
            std::swap(m[i], m[pivot_row]);
            swaps++;
        }

        for (int j = i + 1; j < N; ++j) {
            T factor = m[j][i] / m[i][i];
            for (int k = i; k < N; ++k) {
                m[j][k] -= factor * m[i][k];
            }
        }
    }

    T determinant = 1;
    for (int i = 0; i < N; ++i) {
        determinant *= m[i][i];
    }

    return (swaps % 2 == 0) ? determinant : -determinant;
}

template<int N, typename T> inline T det(mat<N,N,T> m) {
    if (N < 6)  return detRecursive(m);
    else        return detGauss(m);
}

// matrix define
template<int nrows, int ncols, typename T> struct mat {
    vec<ncols,T> rows[nrows] = {{}};

    vec<ncols,T>& operator[] (const int idx)       { assert(idx >= 0 && idx < nrows); return rows[idx]; }
    const vec<ncols,T>& operator[] (const int idx) const { assert(idx >= 0 && idx < nrows); return rows[idx]; }

    // get zero matrix
    static mat<nrows,ncols,T> zero() {
        mat<nrows,ncols,T> res;
        for (int i = nrows; i--; )
            for (int j = ncols; j--; res[i][j] = 0);
        return res;
    }

    // get identity matrix
    static mat<nrows,ncols,T> identity() {
        mat<nrows,ncols,T> res;
        for (int i = nrows; i--; )
            for (int j = ncols; j--; res[i][j] = (i == j));
        return res;
    }

    inline vec<nrows,T> col(const int idx) const {
        assert(idx >= 0 && idx < ncols);
        vec<nrows,T> res;
        for (int i = nrows; i--; res[i] = rows[i][idx]);
        return res;
    }

    inline void setCol(const int idx, const vec<nrows,T>& v) {
        assert(idx >= 0 && idx < ncols);
        for (int i = nrows; i--; rows[i][idx] = v[i]);
    }

    // transpose the matrix
    inline mat<ncols,nrows,T> transpose() const {
        mat<ncols,nrows,T> res;
        for (int i = ncols; i--; res[i] = this->col(i));
        return res;
    }

    // get the complementary submatrix
    inline mat<nrows-1,ncols-1,T> getMinor(const int row, const int col) const {
        mat<nrows-1,ncols-1,T> res;
        for (int i = nrows - 1; i-- ;)
            for (int j = ncols - 1; j--; res[i][j] = rows[i < row ? i : i+1][j < col ? j : j+1]);
        return res;
    }

    // compute the cofactor
    inline T cofactor(const int row, const int col) const {
        return det(getMinor(row, col)) * ((row + col)%2 ? -1 : 1);
    }

    // compute the adjugate matrix
    inline mat<nrows,ncols,T> adjugate() const {
        mat<nrows,ncols,T> res;
        for (int i = nrows; i--; )
            for (int j = ncols; j--; res[j][i] = cofactor(i, j));
        return res;
    }

    // compute the transposed invert matrix
    inline mat<nrows,ncols,T> invertTranspose() const {
        return invert().transpose();
    }

    // compute the invert matrix using Gauss-Jordan elimination
    mat<nrows,ncols,T> invert() const {
        if (nrows != ncols) return zero();
        mat<nrows, nrows*2, T> augmentedMatrix;

        // Construct an augmented matrix [matrix | I]
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nrows; ++j) {
                augmentedMatrix[i][j] = rows[i][j];
                augmentedMatrix[i][j + nrows] = (i == j) ? 1 : 0;
            }
        }

        // Use Gauss-Jordan elimination to change the left part into an identity matrix
        for (int i = 0; i < nrows; ++i) {
            // Set the diagonal element to 1
            T pivot = augmentedMatrix[i][i];
            for (int j = 0; j < nrows * 2; ++j)
                augmentedMatrix[i][j] /= pivot;

            // elimination
            for (int k = 0; k < nrows; ++k) {
                if (k != i) {
                    T factor = augmentedMatrix[k][i];
                    for (int j = 0; j < nrows * 2; ++j)
                        augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                }
            }
        }

        // Extract the right part as the inverse matrix
        mat<nrows, nrows, T> inverseMatrix;
        for (int i = 0; i < nrows; ++i)
            for (int j = 0; j < nrows; ++j)
                inverseMatrix[i][j] = augmentedMatrix[i][j + nrows];

        return inverseMatrix;
    }
};

// override operator <<, print matrix
template<int nrows, int ncols, typename T>
std::ostream& operator<< (std::ostream& out, const mat<nrows,ncols,T>& m) {
    for (int i = 0; i < nrows; i++) std::cout << m[i] << std::endl;
    return out;
}

// override operator +, matrix - matrix addition
template<int nrows, int ncols, typename T> mat<nrows,ncols,T>
inline operator+ (const mat<nrows,ncols,T>& m, const mat<nrows,ncols,T>& n) {
    mat<nrows,ncols,T> res;
    for (int i = nrows; i--; res[i] = m[i] + n[i]) ;
    return res;
}

// override operator +, matrix - scalar addition
template<int nrows, int ncols, typename T> mat<nrows,ncols,T>
inline operator+ (const mat<nrows,ncols,T>& m, const T& k) {
    mat<nrows,ncols,T> res;
    for (int i = nrows; i--; res[i] = m[i] + k);
    return res;
}

// override operator -, matrix - matrix subtraction
template<int nrows, int ncols, typename T> mat<nrows,ncols,T>
inline operator- (const mat<nrows,ncols,T>& m, const mat<nrows,ncols,T>& n) {
    mat<nrows,ncols,T> res;
    for (int i = nrows; i--; res[i] = m[i] - n[i]);
    return res;
}

// override operator -, matrix - scalar addition
template<int nrows, int ncols, typename T> mat<nrows,ncols,T>
inline operator- (const mat<nrows,ncols,T>& m, const T& k) {
    mat<nrows,ncols,T> res;
    for (int i = nrows; i--; res[i] = m[i] - k) ;
    return res;
}

// override operator *, matrix - vector multiplication
template<int nrows, int ncols, typename T>
inline vec<nrows,T> operator* (const mat<nrows,ncols,T>& m, const vec<ncols,T>& v) {
    vec<nrows,T> res;
    for (int i = nrows; i-- ; res[i] = m[i] * v) ;
    return res;
}

// override operator *, matrix - matrix multiplication
template<int M, int N, int K, typename T>
inline mat<M,K,T>operator* (const mat<M,N,T>& lhs, const mat<N,K,T>& rhs) {
    mat<M,K,T> res;
    for (int i = M; i--; )
        for (int j = K; j--; res[i][j] = lhs[i] * rhs.col(j)) ;
    return res;
}

// override operator *, matrix - scalar multiplicat number
template<int nrows, int ncols, typename T>
inline mat<nrows,ncols,T> operator* (const T& k, const mat<nrows,ncols,T>& m) {
    mat<nrows,ncols,T> res;
    for (int i = nrows; i--; res[i] = m[i] * k) ;
    return res;
}

// override operator *, matrix - scalar multiplicat number
template<int nrows, int ncols, typename T>
inline mat<nrows,ncols,T> operator* (const mat<nrows,ncols,T>& m, const T& k) {
    return k * m;
}

// override operator /, matrix - scalar division
template<int nrows, int ncols, typename T>
inline mat<nrows,ncols,T> operator/ (const mat<nrows,ncols,T>& m, const T& k) {
    mat<nrows,ncols,T> res;
    for (int i = nrows; i--; res[i] = m[i] / k) ;
    return res;
}

// Specialization for nrows == ncols == 3
template<typename T> struct mat<3,3,T> {
    vec<3,T> rows[3] = {{}};

    vec<3,T>& operator[] (const int idx)       { assert(idx >= 0 && idx < 3); return rows[idx]; }
    const vec<3,T>& operator[] (const int idx) const { assert(idx >= 0 && idx < 3); return rows[idx]; }

    inline vec<3,T> col(const int idx) const {
        assert(idx >= 0 && idx < 3);
        return { rows[0][idx], rows[1][idx], rows[2][idx] };
    }

    inline void setCol(const int idx, const vec<3,T>& v) {
        assert(idx >= 0 && idx < 3);
        rows[0][idx] = v[0];
        rows[1][idx] = v[1];
        rows[2][idx] = v[2];
    }

    inline mat<3,3,T> transpose() const {
        mat<3, 3, T> res;
        res[0] = this->col(0);
        res[1] = this->col(1);
        res[2] = this->col(2);
        return res;
    }
};

// Specialization for nrows == ncols == 4
template<typename T> struct mat<4,4,T> {
    vec<4,T> rows[4] = {{}};

    vec<4,T>& operator[] (const int idx)       { assert(idx >= 0 && idx < 4); return rows[idx]; }
    const vec<4,T>& operator[] (const int idx) const { assert(idx >= 0 && idx < 4); return rows[idx]; }

    inline vec<4,T> col(const int idx) const {
        assert(idx >= 0 && idx < 4);
        return { rows[0][idx], rows[1][idx], rows[2][idx], rows[3][idx] };
    }

    inline void setCol(const int idx, const vec<4,T>& v) {
        assert(idx >= 0 && idx < 4);
        rows[0][idx] = v[0];
        rows[1][idx] = v[1];
        rows[2][idx] = v[2];
        rows[3][idx] = v[3];
    }

    inline mat<4,4,T> transpose() const {
        mat<4,4,T> res;
        res[0] = this->col(0);
        res[1] = this->col(1);
        res[2] = this->col(2);
        res[3] = this->col(3);
        return res;
    }
};

typedef mat<3,3,double> mat3;
typedef mat<4,4,double> mat4;

inline vec3 operator* (const mat3& m, const vec3& v) {
    return {m[0]*v, m[1]*v, m[2]*v};
}

inline vec4 operator* (const mat4& m, const vec4& v) {
    return {m[0]*v, m[1]*v, m[2]*v, m[3]*v};
}

// Const value defination
static const vec3 XA = {1, 0, 0};
static const vec3 YA = {0, 1, 0};
static const vec3 ZA = {0, 0, 1};
static const vec3 ONE_VEC3 = {1, 1, 1};
static const vec3 ZERO_VEC3 = {0, 0, 0};
static const mat3 E33 = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
static const mat4 E44 = {{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}};

/********************************************************************************
*                                Utility Methods                                *
********************************************************************************/
// Convert degree to radian
inline double deg2rad(double degrees) {
    return degrees * PI / 180.0;
}

// Clamp x to [min, max], default to [0, 1]
inline double clamp(double x, double min=0, double max=1) {
    return (x < min) ? min : ((x > max) ? max : x);
}

// Clamp vec3 to [min, max], default to [0, 1]
inline vec3 clampVec3(const vec3& v, double min=0, double max=1) {
    return {clamp(v.x, min, max), clamp(v.y, min, max), clamp(v.z, min, max)};
}

// Gamma correction with gamma level
inline void gammaCorrection(vec3 &color, double gamma = 2.2) {
    double gamma_recip = 1.0 / gamma;
    color.x = pow(color.x, gamma_recip);
    color.y = pow(color.y, gamma_recip);
    color.z = pow(color.z, gamma_recip);
}

// Generate a random real in [0,1].
inline double randDouble() {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

// Generate a random real in [min,max].
inline double randDouble(double min, double max) {
    return min + (max-min)*randDouble();
}

// Sample in [0, 1]^2
static vec2 sampleUnitSpuare() {
    double x = randDouble();
    double y = randDouble();
    return vec2(x, y);
}

// Sample in uniu disk
static vec2 sampleUnitDisk() {
    double r = randDouble();
    double phi = 2*PI*randDouble();
    return vec2(r*cos(phi), r*sin(phi));
}

static vec3 randVec3() {
    return vec3(randDouble(), randDouble(), randDouble());
}

static vec3 randVec3(double min, double max) {
    return vec3(randDouble(min, max), randDouble(min, max), randDouble(min, max));
}

// Generate random vector in sphere
static vec3 randVecSphere() {
    // Assume that n is normalized
    // double phi = 2*PI*randDouble(), r = randDouble(-1, 1);
    // double a, b;
    // if (r > 0) { a = sqrt(1 - r), b = sqrt(r); }
    // else {a = sqrt(1 + r), b = sqrt(-r); }
    // return vec3(cos(phi)*a, sin(phi)*a, b).normalize();
    double phi = 2*PI*randDouble(), theta = PI*randDouble(-1, 1);
    double a = cos(theta), b = sin(theta);
    return vec3(cos(phi)*a, sin(phi)*a, b).normalize();
}

// Generate random vector in hemisphere
static vec3 randVecSemisphere(const vec3& n) {
    // Assume that n is normalized
    double phi = 2*PI*randDouble(), r = randDouble();
    double a = sqrt(r), b = sqrt(1 - r);
    vec3 u = YA.cross(n);
    u = (u.norm() > EPS ? u : XA.cross(n)).normalize();
    // u = (u.norm() > EPS ? u : XA.cross(n));
    vec3 v = n.cross(u);
    return (u*cos(phi)*a + v*sin(phi)*a + n*b).normalize();
}

// Compute Reflect Ray
inline vec3 reflect(const vec3& v, const vec3& n) {
    return (v - 2*(v*n)*n).normalize();
}

// Compute Refract Ray
inline vec3 refract(const vec3& rin, const vec3& n, double etai_over_etat) {
    double cos_theta1 = rin * n;
    double cos2_theta2 = 1 - etai_over_etat*etai_over_etat*(1 - cos_theta1*cos_theta1);
    if (cos2_theta2 < 0) return reflect(rin, n);    // Total reflection
    vec3 rout_perp =  etai_over_etat * (rin - cos_theta1*n);
    vec3 rout_para = -n * sqrt(cos2_theta2);
    return (rout_perp + rout_para).normalize();
}

inline mat4 scale(double s) {
    mat4 S = E44;
    for (int i = 0; i < 3; i ++) S[i][i] = s;
    return S;
}

inline mat4 translate(const vec3 &t) {
    mat4 T = E44;
    for (int i = 0; i < 3; i ++) T[i][3] = t[i];
    return T;
}

inline mat4 rotateX(double degree) {
    double radian = deg2rad(degree);
    mat4 R;
    R = {{{1.f, 0.f, 0.f, 0.f},
                {0.f, cos(radian), -sin(radian), 0.f},
                {0.f, sin(radian), cos(radian), 0.f},
                {0.f, 0.f, 0.f, 1.f}}};
    return R;
}

inline mat4 rotateY(double degree) {
    double radian = deg2rad(degree);
    mat4 R;
    R = {{{cos(radian), 0.f, sin(radian), 0.f},
                {0.f, 1.f, 0.f, 0.f},
                {-sin(radian), 0.f, cos(radian), 0.f},
                {0.f, 0.f, 0.f, 1.f}}};
    return R;
}

inline mat4 rotateZ(double degree) {
    double radian = deg2rad(degree);
    mat4 R;
    R = {{{ cos(radian), -sin(radian), 0.f, 0.f},
                { sin(radian),  cos(radian), 0.f, 0.f},
                { 0.f, 0.f, 1.f, 0.f},
                { 0.f, 0.f, 0.f, 1.f}}};
    return R;
}

inline mat4 rotate(const vec3 &n, double degree) {
    mat4 R = E44;
    return R;
}

#endif