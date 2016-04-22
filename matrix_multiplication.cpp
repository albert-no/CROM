/* 
    matrix_multiplication.cpp
    butterfly matrix multiplication
*/
#include "matrix_multiplication.hpp"

void butterfly_matrix_multiplication (double *x,
                                      double *thetas,
                                      int half_mat_dim,
                                      int x_start_idx,
                                      int x_end_idx,
                                      int theta_start_idx,
                                      int theta_end_idx,
                                      int mat_idx) {
    /*
       Compute butterfly matrix - vector multiplication
       Recursively compute based on mat_idx
       dimension should be the power of 2

       Input Parameters
       ----------------
       x :: input vector
       thetas :: array of thetas
       x_idx :: start and end indices of x
       theta_idx :: start and end indices of theta
       mat_idx :: category of matrix index
    */

    double old_x[2];
    double new_x[2];
    double theta, costheta, sintheta;

    int iter_idx;
    int x_mid_idx = (x_start_idx+x_end_idx)/2;
    int theta_mid_idx;
    // mat_idx = 0 
    if (mat_idx == 0) {
        for (iter_idx=x_start_idx; iter_idx<x_mid_idx; iter_idx++) {
            old_x[0] = x[iter_idx];
            old_x[1] = x[iter_idx+half_mat_dim];
            theta = thetas[iter_idx];
            costheta = cos(theta);
            sintheta = sin(theta);
            new_x[0] = costheta * old_x[0] - sintheta * old_x[1];
            new_x[1] = sintheta * old_x[0] + costheta * old_x[1];
            x[iter_idx] = new_x[0];
            x[iter_idx+half_mat_dim] = new_x[1];
        }
    }
    // mat_idx > 0
    else {
        theta_mid_idx = (theta_start_idx+theta_end_idx)/2;
        half_mat_dim /= 2;
        butterfly_matrix_multiplication(x,
                                        thetas,
                                        half_mat_dim,
                                        x_start_idx,
                                        x_mid_idx,
                                        theta_start_idx,
                                        theta_mid_idx,
                                        mat_idx-1);
        butterfly_matrix_multiplication(x,
                                        thetas,
                                        half_mat_dim,
                                        x_mid_idx,
                                        x_end_idx,
                                        theta_mid_idx,
                                        theta_end_idx,
                                        mat_idx-1);
    }
}
