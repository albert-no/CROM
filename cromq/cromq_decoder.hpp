// cromq_decoder.hpp
#ifndef CROMQ_DECODER_H
#define CROMQ_DECODER_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/SVD>

#include "cromq_util.hpp"
#include "../decoder/crom_decoder.hpp"

/* 
*/

class CROMq_decoder
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
    double R_dec;
    bool verbose;

    // Extract from original q scores
    std::vector<double> mu;

    // Result of SVD of covariance matrix
    std::vector<std::vector<double>> v_mat;
    std::vector<double> std_array;

    // x_hat before unnormalize
    std::vector<std::vector<double>> x_hat_array;

    // reconstructed q_scores
    std::vector<std::vector<char>> q_scores;

    // Allocated rates
    std::vector<double> R_array;

    // round q_score into char format
    char q_score_round(double q_score);

    // normalize q_scores with mu vector and V matrix
    void unnormalize_q_scores();

    // write q_scores
    void write_q_scores();


public:
    // Constructor
    CROMq_decoder(std::string name_in,
                  std::string fname_in,
                  int id_in,
                  int num_x_in,
                  int x_dim_in,
                  double rd_param_in,
                  double R_dec_in,
                  std::vector<double> R_allocated_in,
                  std::vector<double> mu_in,
                  std::vector<double> std_array_in,
                  std::vector<std::vector<double>> v_mat_in,
                  bool verbose);

    // Destructor
    ~CROMq_decoder();

    // Run CROMq_decoder
    void run();
};

#endif
