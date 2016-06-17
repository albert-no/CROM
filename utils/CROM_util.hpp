/*
   CROM_util.hpp
*/
#ifndef CROM_UTIL_H
#define CROM_UTIL_H

#include <cmath>
#include <iostream>
#include <stdio.h>

int find_max_index(double *x, int x_dim);
double compute_l2(double *x, int x_dim);
double compute_l2_dist(double *x, double *x_hat, int x_dim);
void print_vector(double *x, int x_dim);
void unnormalize_vector(double *x, int x_dim);
void copy_vector(double *x, double *xout, int x_dim);
void normalize_then_copy_vector(double *x, double *xout, int x_dim);

#endif
