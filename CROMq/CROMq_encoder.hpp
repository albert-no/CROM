// CROMq_encoder.hpp
#ifndef CROMQ_ENCODER_H
#define CROMQ_ENCODER_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "../encoder/CROM_encoder.hpp"

/* 
*/

class CROMq_encoder
{
    /* read file and split rate and run CROM

    Parameters
    ----------
    name :: name of the object
    fname :: name of the file contains sequences
    num_x :: number of sequences, default = 101
    R_array :: num_x dimension array that contains rate of subsequences

    */
    std::string name;
    std::string fname;
    int num_x;
    int x_dim;
    double *std_array;
    double *R_array;



public:
    // Constructor
    CROMq_encoder(std::string name_in, std::string fname_in, int num_x_in,
                  int x_dim_in, double *std_array_in, double *R_array_in);

    // Destructor
    ~CROMq_encoder();

    // Run CROMq_encoder
    void run();
}

#endif
