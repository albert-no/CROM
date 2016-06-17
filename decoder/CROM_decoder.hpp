/* 
   CROM_decoder.hpp
*/
#ifndef CROM_DECODER_H
#define CROM_DECODER_H

#include "../utils/matrix_multiplication.hpp"
#include "../utils/CROM_util.hpp"
#include <cstdlib>
#include <cmath>

void CROM_decoding_step(double *x, int x_dim, double scale, int m);
void CROM_decoder(double *xhat, int x_dim, int L, int *m_array, bool verbose);

#endif
