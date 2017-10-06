/* 
   CROM_decoder.hpp
*/
#ifndef CROM_DECODER_H
#define CROM_DECODER_H

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>

#include "../utils/matrix_multiplication.hpp"
#include "../utils/crom_util.hpp"

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
        verbose :: whether printing intermediate values
        name :: name of the object
    */

    double *x_hat;
    int x_dim;
    int L;
    int *m_array;
    bool verbose;
    std::string name;

    // Single iteration of CROM_decoder with k=1
    // scale :: scale factor of iteration
    // m :: message (index of maximum element) of iteration
    void step(double scale, int m);

public:
    // Constructor
    CROM_decoder(std::string name_in, int x_dim_in, int L_in, bool verbose_in);

    // Destructor
    ~CROM_decoder();

    // set the array of messages (m_array) via copying from m_array_in
    void set_m_array(int *m_array_in);

    // read the array of messages (m_array) from either binary or txt file
    void read_m_array(bool binary);

    // run CROM decoder
    void run();

    // read x_hat via copying to x_hat_copy
    void copy_x_hat(double *x_hat_copy);

    // write x_hat
    void write_x_hat();
};

#endif
