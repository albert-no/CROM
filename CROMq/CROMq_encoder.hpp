// CROMq_encoder.hpp
#ifndef CROMQ_ENCODER_H
#define CROMQ_ENCODER_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include <Eigen/SVD>

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
    x_dim :: dimension of each sequence

    std_array :: num_x dimension array that contains std of subsequences
    R_array :: num_x dimension array that contains rate of subsequences
    x_array :: x_dim * num_x dimension array that contains subsequences
    */

    std::string name;
    std::string fname;
    int num_x;
    int x_dim;

    // Extract from original q scores
    double* mu;
    double** cov;

    // Result of SVD of covariance matrix
    double** v_mat;
    double* std_array;

    // Normalized q scores
    double** x_array;

    // Allocated rates
    double* R_array;


    // Allocate rate according to the std_array
    void allocate_rate();

    // Get x_array();
    void get_x_array();

    // Get covariance matrix
    void get_cov();


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
