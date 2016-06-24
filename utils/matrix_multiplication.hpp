/* 
    matrix_multiplication.hpp
    butterfly matrix multiplication
*/
#ifndef MATRIX_MULTIPLICATION_H
#define MATRIX_MULTIPLICATION_H

#include <cmath>
#include <fftw3.h>

void butterfly_matrix_multiplication (double *x,
                                      double *thetas,
                                      int half_len,
                                      int x_start_idx,
                                      int theta_start_idx,
                                      int mat_idx);
#endif
