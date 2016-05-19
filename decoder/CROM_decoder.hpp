/* 
   CROM_decoder.hpp
*/
#ifndef CROM_DECODER_H
#define CROM_DECODER_H

#include "../utils/matrix_multiplication.hpp"
#include "../utils/CROM_util.hpp"
#include <cstdlib>
#include <cmath>

#define THETA_SEED 1234

int CROM_decoding_step(double *x, int x_dim, double scale);
void CROM_decoder(double *x, int x_dim, int L, int *m_array, bool print_l2norm);

#endif
