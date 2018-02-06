/* 
   CROM_decoder.hpp
*/
#ifndef CROM_DECODER_H
#define CROM_DECODER_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

#include "../utils/CROM_util.hpp"
#include "../utils/matrix_multiplication.hpp"


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

    int L;
    int x_dim;
    bool verbose;

    std::string name;
    std::vector<double> x_hat;
    std::vector<int> m_array;

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
    void set_m_array(std::vector<int> &m_array_in);

    // read the array of messages (m_array) from either binary or txt file
    void read_m_array(bool binary);

    // run CROM decoder
    void run();

    // read x_hat via copying to x_hat_copy
    void copy_x_hat(std::vector<double> &x_hat_copy);

    // write x_hat
    void write_x_hat();
};

#endif
