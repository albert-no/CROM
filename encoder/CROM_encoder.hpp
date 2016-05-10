/* 
   CROM_encoder.hpp
*/
#ifndef CROM_ENCODER_H
#define CROM_ENCODER_H

#include "../utils/matrix_multiplication.hpp"
#include <cstdlib>
#include <cmath>

#define THETA_SEED 1234

int find_max_index(double *x, int x_dim);
int CROM_step(double *x, int x_dim, double scale);
void CROM_encoder(double *x, int x_dim, int L, int *m_array);

#endif
