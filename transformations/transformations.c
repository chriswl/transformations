/* transformations.c

A Python C extension module for homogeneous transformation matrices and
quaternions.

Refer to the transformations.py module for documentation and tests.

:Author:
  `Christoph Gohlke <http://www.lfd.uci.edu/~gohlke/>`_

:Organization:
  Laboratory for Fluorescence Dynamics, University of California, Irvine

:Version: 2015.07.18

Install
-------
Use this Python distutils setup script to build the extension module::

  # setup.py
  # Usage: ``python setup.py build_ext --inplace``
  from distutils.core import setup, Extension
  import numpy
  setup(name='_transformations',
        ext_modules=[Extension('_transformations', ['transformations.c'],
                               include_dirs=[numpy.get_include()])])

License
-------
Copyright (c) 2007-2015, Christoph Gohlke
Copyright (c) 2007-2015, The Regents of the University of California
Produced at the Laboratory for Fluorescence Dynamics
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
* Neither the name of the copyright holders nor the names of any
  contributors may be used to endorse or promote products derived
  from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#define _VERSION_ "2015.07.18"

#define WIN32_LEAN_AND_MEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Python.h"

#ifdef _WIN32
#include <windows.h>
#include <wincrypt.h>
#endif

#include "structmember.h"
#include "math.h"
#include "float.h"
#include "string.h"
#include "numpy/arrayobject.h"

#define EPSILON 8.8817841970012523e-016 /* 4.0 * DBL_EPSILON */
#define PIVOT_TOLERANCE 1.0e-14
#define DEG2RAD 0.017453292519943295769236907684886
#define TWOPI 6.283185307179586476925286766559
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

#define ISZERO(x) (((x) < EPSILON) && ((x) > -EPSILON))
#define NOTZERO(x) (((x) > EPSILON) || ((x) < -EPSILON))

/*****************************************************************************/
/* C helper functions */

/*
Return random doubles in half-open interval [0.0, 1.0).
Uses /dev/urandom or CryptoAPI. Not very fast but random.
Assumes sizeof(double) == 2*sizeof(int).
*/
int random_doubles(
    double *buffer,
    Py_ssize_t size)
{
#ifndef _WIN32
    int done;
    FILE *rfile;
    if (size < 1)
        return 0;
    rfile = fopen("/dev/urandom", "rb");
    if (rfile == NULL)
        return -1;
    done = fread(buffer, size*sizeof(double), 1, rfile);
    fclose(rfile);
#else
#pragma comment(lib,"advapi32")
    BOOL done;
    HCRYPTPROV hcryptprov;
    if (size < 1)
        return 0;
    if (!CryptAcquireContext(&hcryptprov, NULL, NULL, PROV_RSA_FULL,
                             CRYPT_VERIFYCONTEXT) || !hcryptprov)
        return -1;
    done = CryptGenRandom(hcryptprov, (DWORD)(size*sizeof(double)),
                          (unsigned char *)buffer);
    CryptReleaseContext(hcryptprov, 0);
#endif
    if (done) {
        unsigned int a, b;
        unsigned int *p = (unsigned int *)buffer;
        while (size--) {
            a = (*p++) >> 5;
            b = (*p++) >> 6;
            *buffer++ = (a * 67108864.0 + b) / 9007199254740992.0;
        }
        return 0;
    }
    return -1;
}

/*
Tridiagonal matrix from symmetric 4x4 matrix using Housholder reduction.
The input matrix is altered.
*/
int tridiagonalize_symmetric_44(
    double *matrix,      /* double[16] */
    double *diagonal,    /* double[4] */
    double *subdiagonal) /* double[3] */
{
    double t, n, g, h, v0, v1, v2;
    double *u;
    double *M = matrix;

    u = &M[1];
    t = u[1]*u[1] + u[2]*u[2];
    n = sqrt(u[0]*u[0] + t);
    if (n > EPSILON) {
        if (u[0] < 0.0)
            n = -n;
        u[0] += n;
        h = (u[0]*u[0] + t) / 2.0;
        v0 = M[5]*u[0] + M[6]*u[1]  + M[7]*u[2];
        v1 = M[6]*u[0] + M[10]*u[1] + M[11]*u[2];
        v2 = M[7]*u[0] + M[11]*u[1] + M[15]*u[2];
        v0 /= h;
        v1 /= h;
        v2 /= h;
        g = (u[0]*v0 + u[1]*v1 + u[2]*v2) / (2.0 * h);
        v0 -= g * u[0];
        v1 -= g * u[1];
        v2 -= g * u[2];
        M[5] -=  2.0*v0*u[0];
        M[10] -= 2.0*v1*u[1];
        M[15] -= 2.0*v2*u[2];
        M[6]  -= v1*u[0] + v0*u[1];
        M[7]  -= v2*u[0] + v0*u[2];
        M[11] -= v2*u[1] + v1*u[2];
        M[1] = -n;
    }

    u = &M[6];
    t = u[1]*u[1];
    n = sqrt(u[0]*u[0] + t);
    if (n > EPSILON) {
        if (u[0] < 0.0)
            n = -n;
        u[0] += n;
        h = (u[0]*u[0] + t) / 2.0;
        v0 = M[10]*u[0] + M[11]*u[1];
        v1 = M[11]*u[0] + M[15]*u[1];
        v0 /= h;
        v1 /= h;
        g = (u[0]*v0 + u[1]*v1) / (2.0 * h);
        v0 -= g * u[0];
        v1 -= g * u[1];
        M[10] -= 2.0*v0*u[0];
        M[15] -= 2.0*v1*u[1];
        M[11] -= v1*u[0] + v0*u[1];
        M[6] = -n;
    }

    diagonal[0] = M[0];
    diagonal[1] = M[5];
    diagonal[2] = M[10];
    diagonal[3] = M[15];
    subdiagonal[0] = M[1];
    subdiagonal[1] = M[6];
    subdiagonal[2] = M[11];

    return 0;
}

/*
Return largest eigenvalue of symmetric tridiagonal matrix.
Matrix Algorithms: Volume II: Eigensystems. By GW Stewart. 
Chapter 3. page 197.
*/
double max_eigenvalue_of_tridiag_44(
    double *diagonal,    /* double[4] */
    double *subdiagonal) /* double[3] */
{
    int count;
    double lower, upper, t0, t1, d, eps, eigenv;
    double *a = diagonal;
    double *b = subdiagonal;

    /* upper and lower bounds using Gerschgorin's theorem */
    t0 = fabs(b[0]);
    t1 = fabs(b[1]);
    lower = a[0] - t0;
    upper = a[0] + t0;
    d = a[1] - t0 - t1;
    lower = MIN(lower, d);
    d = a[1] + t0 + t1;
    upper = MAX(upper, d);
    t0 = fabs(b[2]);
    d = a[2] - t0 - t1;
    lower = MIN(lower, d);
    d = a[2] + t0 + t1;
    upper = MAX(upper, d);
    d = a[3] - t0;
    lower = MIN(lower, d);
    d = a[3] + t0;
    upper = MAX(upper, d);

    /* precision */
    /* eps = (4.0 * (fabs(lower) + fabs(upper))) * DBL_EPSILON; */
    eps = 1e-18;

    /* interval bisection until width is less than tolerance */
    while (fabs(upper - lower) > eps) {

        eigenv = (upper + lower) / 2.0;

        if ((eigenv == upper) || (eigenv == lower))
            return eigenv;

        /* counting pivots < 0 */
        d = a[0] - eigenv;
        count = (d < 0.0) ? 1 : 0;
        if (fabs(d) < eps)
            d = eps;
        d = a[1] - eigenv - b[0]*b[0] / d;
        if (d < 0.0)
            count++;
        if (fabs(d) < eps)
            d = eps;
        d = a[2] - eigenv - b[1]*b[1] / d;
        if (d < 0.0)
            count++;
        if (fabs(d) < eps)
            d = eps;
        d = a[3] - eigenv - b[2]*b[2] / d;
        if (d < 0.0)
            count++;

        if (count < 4)
            lower = eigenv;
        else
            upper = eigenv;
    }

    return (upper + lower) / 2.0;
}

/*
Eigenvector of symmetric tridiagonal 4x4 matrix using Cramer's rule.
*/
int eigenvector_of_symmetric_44(
    double *matrix, /* double[16] */
    double *vector, /* double[4]  */
    double *buffer) /* double[12] */
{
    double n, eps;
    double *M = matrix;
    double *v = vector;
    double *t = buffer;

    /* eps: minimum length of eigenvector to use */
    eps = (M[0]*M[5]*M[10]*M[15] - M[1]*M[1]*M[11]*M[11]) * 1e-6;
    eps *= eps;
    if (eps < EPSILON)
        eps = EPSILON;

    t[0] = M[10] * M[15];
    t[1] = M[11] * M[11];
    t[2] = M[6] *  M[15];
    t[3] = M[11] * M[7];
    t[4] = M[6] *  M[11];
    t[5] = M[10] * M[7];
    t[6] = M[2] *  M[15];
    t[7] = M[11] * M[3];
    t[8] = M[2] *  M[11];
    t[9] = M[10] * M[3];
    t[10] = M[2] * M[7];
    t[11] = M[6] * M[3];

    v[0] =  t[1]*M[1] + t[6]*M[6] + t[9]*M[7];
    v[0] -= t[0]*M[1] + t[7]*M[6] + t[8]*M[7];
    v[1] =  t[2]*M[1] + t[7]*M[5] + t[10]*M[7];
    v[1] -= t[3]*M[1] + t[6]*M[5] + t[11]*M[7];
    v[2] =  t[5]*M[1] + t[8]*M[5] + t[11]*M[6];
    v[2] -= t[4]*M[1] + t[9]*M[5] + t[10]*M[6];
    v[3] =  t[0]*M[5] + t[3]*M[6] + t[4]*M[7];
    v[3] -= t[1]*M[5] + t[2]*M[6] + t[5]*M[7];
    n = v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3];

    if (n < eps) {
        v[0] =  t[0]*M[0] + t[7]*M[2] + t[8]*M[3];
        v[0] -= t[1]*M[0] + t[6]*M[2] + t[9]*M[3];
        v[1] =  t[3]*M[0] + t[6]*M[1] + t[11]*M[3];
        v[1] -= t[2]*M[0] + t[7]*M[1] + t[10]*M[3];
        v[2] =  t[4]*M[0] + t[9]*M[1] + t[10]*M[2];
        v[2] -= t[5]*M[0] + t[8]*M[1] + t[11]*M[2];
        v[3] =  t[1]*M[1] + t[2]*M[2] + t[5]*M[3];
        v[3] -= t[0]*M[1] + t[3]*M[2] + t[4]*M[3];
        n = v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3];
    }

    if (n < eps) {
        t[0]  = M[2] * M[7];
        t[1]  = M[3] * M[6];
        t[2]  = M[1] * M[7];
        t[3]  = M[3] * M[5];
        t[4]  = M[1] * M[6];
        t[5]  = M[2] * M[5];
        t[6]  = M[0] * M[7];
        t[7]  = M[3] * M[1];
        t[8]  = M[0] * M[6];
        t[9]  = M[2] * M[1];
        t[10] = M[0] * M[5];
        t[11] = M[1] * M[1];

        v[0] =  t[1]*M[3] + t[6]*M[11] + t[9]*M[15];
        v[0] -= t[0]*M[3] + t[7]*M[11] + t[8]*M[15];
        v[1] =  t[2]*M[3] + t[7]*M[7] + t[10]*M[15];
        v[1] -= t[3]*M[3] + t[6]*M[7] + t[11]*M[15];
        v[2] =  t[5]*M[3] + t[8]*M[7] + t[11]*M[11];
        v[2] -= t[4]*M[3] + t[9]*M[7] + t[10]*M[11];
        v[3] =  t[0]*M[7] + t[3]*M[11] + t[4]*M[15];
        v[3] -= t[1]*M[7] + t[2]*M[11] + t[5]*M[15];
        n = v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3];
    }

    if (n < eps) {
        v[0] =  t[8]*M[11] + t[0]*M[2] + t[7]*M[10];
        v[0] -= t[6]*M[10] + t[9]*M[11] + t[1]*M[2];
        v[1] =  t[6]*M[6] + t[11]*M[11] + t[3]*M[2];
        v[1] -= t[10]*M[11] + t[2]*M[2] + t[7]*M[6];
        v[2] =  t[10]*M[10] + t[4]*M[2] + t[9]*M[6];
        v[2] -= t[8]*M[6] + t[11]*M[10] + t[5]*M[2];
        v[3] =  t[2]*M[10] + t[5]*M[11] + t[1]*M[6];
        v[3] -= t[4]*M[11] + t[0]*M[6] + t[3]*M[10];
        n = v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3];
    }

    if (n < eps)
        return -1;

    n = sqrt(n);
    v[0] /= n;
    v[1] /= n;
    v[2] /= n;
    v[3] /= n;

    return 0;
}

/*
Matrix 2x2 inversion.
*/
int invert_matrix22(
    double *matrix, /* double[4] */
    double *result) /* double[4] */
{
    double det = matrix[0]*matrix[3] - matrix[1]*matrix[2];

    if (ISZERO(det))
        return -1;

    result[0] = matrix[3] / det;
    result[1] = -matrix[1] / det;
    result[2] = -matrix[2] / det;
    result[3] = matrix[0] / det;

    return 0;
}

/*
Matrix 3x3 inversion.
*/
int invert_matrix33(
    double *matrix, /* double[9] */
    double *result) /* double[9] */
{
    int i;
    double det;
    double *M = matrix;

    result[0] = M[8]*M[4] - M[7]*M[5];
    result[1] = M[7]*M[2] - M[8]*M[1];
    result[2] = M[5]*M[1] - M[4]*M[2];
    result[3] = M[6]*M[5] - M[8]*M[3];
    result[4] = M[8]*M[0] - M[6]*M[2];
    result[5] = M[3]*M[2] - M[5]*M[0];
    result[6] = M[7]*M[3] - M[6]*M[4];
    result[7] = M[6]*M[1] - M[7]*M[0];
    result[8] = M[4]*M[0] - M[3]*M[1];

    det = M[0]*result[0] + M[3]*result[1] + M[6]*result[2];

    if (ISZERO(det))
        return -1;

    det = 1.0 / det;
    for (i = 0; i < 9; i++)
        result[i] *= det;

    return 0;
}

/*
Matrix 4x4 inversion.
*/
int invert_matrix44(
    double *matrix, /* double[16] */
    double *result) /* double[16] */
{
    int i;
    double det;
    double t[12];
    double *M = matrix;

    t[0] = M[10] * M[15];
    t[1] = M[14] * M[11];
    t[2] = M[6] * M[15];
    t[3] = M[14] * M[7];
    t[4] = M[6] * M[11];
    t[5] = M[10] * M[7];
    t[6] = M[2] * M[15];
    t[7] = M[14] * M[3];
    t[8] = M[2] * M[11];
    t[9] = M[10] * M[3];
    t[10] = M[2] * M[7];
    t[11] = M[6] * M[3];

    result[0] = t[0]*M[5] + t[3]*M[9] + t[4]*M[13];
    result[0] -= t[1]*M[5] + t[2]*M[9] + t[5]*M[13];
    result[1] = t[1]*M[1] + t[6]
