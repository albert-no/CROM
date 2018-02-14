// cromq_encoder.hpp
#ifndef CROMQ_ENCODER_H
#define CROMQ_ENCODER_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/SVD>

#include "../encoder/crom_encoder.hpp"
#include "cromq_util.hpp"

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
    rd_param :: assuming D = e^{-rd_param * R}, default = 1.4
    R_overall :: overall rate R of CROMq

    std_array :: num_x dimension array that contains std of subsequences
    R_array :: num_x dimension array that contains rate of subsequences
    x_array :: x_dim * num_x dimension array that contains subsequences
    */

    std::string name;
    std::string fname;
    int id;
    int num_x;
    int x_dim;
    double rd_param;
    double R_overall;
    bool verbose;

    // Extract from original q scores
    std::vector<double> mu;
    std::vector<std::vector<double>> cov;

    // Result of SVD of covariance matrix
    std::vector<std::vector<double>> v_mat;
    std::vector<double> std_array;

    // Normalized q scores
    std::vector<std::vector<double>> x_array;

    // Allocated rates
    std::vector<double> R_array;

    // Number of nonzero rate
    int num_nonzero_rate;

    // read qscores and compute mu
    void get_q_scores_and_mu(std::vector<std::vector<double>> &q_scores);

    // normalize q_scores with mu vector and V matrix
    void normalize_q_scores(std::vector<std::vector<double>> &q_scores);

    // Get x_array();
    void get_x_array();

    // Get covariance matrix
    void get_cov(std::vector<std::vector<double>> &q_scores);

    // Do SVD
    void do_svd();

    // write mu_array, std_array, v_mat
    void write_svd_params();


public:
    // Constructor
    CROMq_encoder(std::string name_in, std::string fname_in, int id_in, int num_x_in,
                  int x_dim_in, double rd_param_in, double R_overall_in, bool verbose);

    // Destructor
    ~CROMq_encoder();

    // Run CROMq_encoder
    void run();
};

#endif
