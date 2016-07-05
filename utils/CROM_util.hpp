// CROM_util.hpp
#ifndef CROM_UTIL_H
#define CROM_UTIL_H

#include <cmath>
#include <iostream>
#include <stdio.h>

/* Find the index of maximum element of the vector x

   Parameters
   ----------
   x :: input vector
   x_dim :: dimension of x

   Returns
   -------
   max_idx :: index of maximum element */
int find_max_index(double *x, int x_dim);

/* Compute l2 norm of the vector x

   Parameters
   ----------
   x :: input vector
   x_dim :: dimension of x

   Returns
   -------
   l2norm :: l2norm of x */
double compute_l2(double *x, int x_dim);

/* Compute l2 distance between two vectors x and x_hat

   Parameters
   ----------
   x :: input vector
   x_hat :: input vector
   x_dim :: dimension of x

   Returns
   -------
   l2norm :: l2 distance of of x and xhat */
double compute_l2_dist(double *x, double *x_hat, int x_dim);

/* Print vector x

   Parameters
   ----------
   x :: input vector
   x_dim :: dimension of x */
void print_vector(double *x, int x_dim);

/* Unnormalize the vector x before idct2

   Reverse of normalize_then_copy_vector

   Parameters
   ----------
   x :: input vector
   x_dim :: dimension of x */
void unnormalize_vector(double *x, int x_dim);

/* Copy vector x to xout

   Parameters
   ----------
   x :: input vector
   x_out :: output vector
   x_dim :: dimension of x */
void copy_vector(double *x, double *xout, int x_dim);

/* Normalize the vector xout and copy to x after dct2

   DCT II of fftw3 is a standard one. In order to normalize it
   (or make DCT-II matrix to be an orthogonal matrix), we need
   to scale X_0 by \sqrt{1/n} and all others by \sqrt{2/N}.

   Parameters
   ----------
   x :: input vector
   x_out :: output vector
   x_dim :: dimension of x */
void normalize_then_copy_vector(double *x, double *xout, int x_dim);

#endif
