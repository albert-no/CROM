/*
   CROM_util.cpp
*/
#include "CROM_util.hpp"

int find_max_index(double *x, int x_dim){
    /*
       Function that returns the index of maximum element

       Input Parameters
       ______________-_
       x :: input vector
       x_dim :: dimension of x
       scale :: scale factor of iteration

       Return Parameters
       _________________
       max_idx :: index of maximum element
    */

    int iter_idx;
    int max_idx = 0;
    double max_val = x[0];

    for (iter_idx=1; iter_idx<x_dim; iter_idx++) {
        if (max_val < x[iter_idx]) {
            max_idx = iter_idx;
            max_val = x[iter_idx];
        }
    }
    return max_idx;
}
