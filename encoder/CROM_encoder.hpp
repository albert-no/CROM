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

//int CROM_step(double *x, int x_dim, double scale);
//void CROM_encoder(double *x, int x_dim, int L, int *m_array, bool verbose);

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

    // encoding step
    int step(double scale);

public:
    CROM_encoder(int x_dim_input, double R_input, bool verbose_input);
    ~CROM_encoder();

    // TBD XXX we may want to update the way of reading x value
    void read_x(std::string filename);
    void run();
    void print_m_array();
    // TBD XXX we may want to print m_array to file
    // need m_array_file_name
    // we may want to have x_file as a class member

};

#endif
