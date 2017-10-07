/* 
   CROM_encoder.hpp
*/
#ifndef CROM_ENCODER_H
#define CROM_ENCODER_H

#include <cmath>
#include <fstream>
#include <iostream>

#include "../utils/crom_util.hpp"
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
        l2_array :: array of l2norms after each step 
        verbose :: whether printing intermediate l2 norm
        name :: name of the object
    */
    double *x;
    int x_dim;
    double R;
    int L;
    int *m_array;
    double *l2_array;
    bool verbose;
    std::string name;

    // Single iteration of CROM_encoder with k=1
    // scale :: scale factor of iteration
    int step(double scale);

public:
    // Constructor
    CROM_encoder(std::string name_in, int x_dim_in, double R_in, bool verbose_in);

    // Destructor
    ~CROM_encoder();

    // set input vector x via copying from x_in
    void set_x(double *x_in);

    // set input vector x via copying from x_array_in
    // this feature is being used in CROMq
    void set_x_from_array(double **x_array_in, int idx, bool vertical);

    // read vector x via copying to x_copy
    void copy_x(double *x_copy);

    // read m_array via copying to m_array_copy
    void copy_m_array(int *m_array_copy);

    // read l2_array via copying to l2_array_copy
    void copy_l2_array(double *l2_array_copy);

    // get number of iterations L
    int get_L();

    // run CROM encoder
    void run();

    // print m_array
    void print_m_array();

    // write m_array
    void write_m_array(bool binary);
};

#endif