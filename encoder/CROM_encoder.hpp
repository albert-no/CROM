/* 
   CROM_encoder.hpp
*/
#ifndef CROM_ENCODER_H
#define CROM_ENCODER_H

#include <cstdlib>
#include <cmath>
#include <fstream>

#include "../utils/CROM_util.hpp"
#include "../utils/matrix_multiplication.hpp"


class CROM_encoder
{
    /* CROM encoder with k=1

        Parameters
        ----------
        x :: input vector
        x_dim :: dimension of x
        R :: rate
        L :: number of iterations
        m_array :: array of massages
        verbose :: whether printing intermediate l2 norm
    */
    double *x;
    double R;
    int x_dim;
    int L;
    int *m_array;
    bool verbose;

    // Single iteration of CROM_encoder with k=1
    // scale :: scale factor of iteration
    int step(double scale);

public:
    // Constructor
    CROM_encoder(int x_dim_input, double R_input, bool verbose_input);

    // Destructor
    ~CROM_encoder();

    // set input vector x via copy from x_in
    void set_x(double *x_in);

    // read vector x via copying to x_copy
    void copy_x(double *x_copy);

    // read m_array via copying to m_array_copy
    void copy_m_array(int *m_array_copy);

    // get number of iterations L
    int get_L();

    // run CROM encoder
    void run();

    // print m_array
    void print_m_array();

    // TBD XXX we may want to print m_array to file
    // need m_array_file_name
    // we may want to have x_file as a class member

};

#endif
