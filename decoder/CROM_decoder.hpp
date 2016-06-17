/* 
   CROM_decoder.hpp
*/
#ifndef CROM_DECODER_H
#define CROM_DECODER_H

#include <cstdlib>
#include <cmath>
#include <fstream>

#include "../utils/matrix_multiplication.hpp"
#include "../utils/CROM_util.hpp"

// void CROM_decoding_step(double *x, int x_dim, double scale, int m);
// void CROM_decoder(double *xhat, int x_dim, int L, int *m_array, bool verbose);

class CROM_decoder {
    /* CROM decoder with k=1 

        Parameters
        ----------
        xhat :: reconstruction vector
            initialized with 0
        x_dim :: dimension of x
        L :: number of iterations
        m_array :: array of massages
    */

    double *x_hat;
    int x_dim;
    int L;
    int *m_array;
    bool verbose;

    void step(double scale, int m);

public:
    CROM_decoder(int x_dim_in, int L_in, bool verbose_in);
    ~CROM_decoder();
    void set_m_array(int *m_array_in);
    void run();
    void copy_x_hat(double *x_hat_copy);

};

#endif
