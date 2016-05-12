/* 
    matrix_multiplication.hpp
    butterfly matrix multiplication
*/
#ifndef MATRIX_MULTIPLICATION_H
#define MATRIX_MULTIPLICATION_H

#include <cmath>
#include <fftw3.h>

#define BLOCKLENGTH 65536
#define HALFBLOCKLENGTH 32768

void butterfly_matrix_multiplication (double *x,
                                      double *thetas,
                                      int half_mat_dim,
                                      int x_start_idx,
                                      int x_end_idx,
                                      int theta_start_idx,
                                      int theta_end_idx,
                                      int mat_idx);
#endif
