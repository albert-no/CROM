// cromq_util.hpp
#ifndef CROMQ_UTIL_H
#define CROMQ_UTIL_H

#include <cmath>
#include <iostream>
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

#endif
