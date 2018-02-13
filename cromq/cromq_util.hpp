// cromq_util.hpp
#ifndef CROMQ_UTIL_H
#define CROMQ_UTIL_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


/* Allocate rate according to the std_array
   Assuming D = e^{-1.4R}

   Parameters
   ----------
   std_array:: array of std
   R_array:: array of rates
   num_x:: number of columns
   rd_param:: controlling parameter of D(R), default is 1.4
   R_overall:: overall target rate

   Returns
   -------
   num_nonzero_rate:: number of nonzero rates in R_array */
int allocate_rate(std::vector<double> &std_array, std::vector<double> &R_array, int num_x,
        double rd_param, double R_overall);

/* Get the file name that has mu values

   Parameters
   ----------
   fname:: string of file name
   id:: id

   Returns
   -------
   mu_fname:: file name that contains mu values */
std::string get_mu_fname(std::string fname, int id);

/* Get the file name that has v_mat values

   Parameters
   ----------
   fname:: string of file name
   id:: id

   Returns
   -------
   v_mat_fname:: file name that contains v_mat values */
std::string get_v_mat_fname(std::string fname, int id);

/* Get the file name that has mu values

   Parameters
   ----------
   fname:: string of file name
   id:: id

   Returns
   -------
   std_array_fname:: file name that contains mu values */
std::string get_std_array_fname(std::string fname, int id);

/* Compute the distortion between original and the reconstructed q_scores

   Parameters
   ----------
   ifname:: file name of original q_score
   ofname:: file name of reconstructed q_score
   num_x:: number of columns
   x_dim:: dimension of x (number of rows)

   Returns
   -------
   distortion:: normalized distortion between original and reconstructed q_scores */
double compute_distortion(std::string ifname, std::string ofname, int num_x, int x_dim);

/* Get the output file name that has (or will have) reconstructed q_scores

   Parameters
   ----------
   name:: name of CROMq decoder
   id:: id
   R_overall:: target rate

   Returns
   -------
   ofname:: output file name */
std::string get_ofname(std::string name, int id, double R_overall);
#endif
