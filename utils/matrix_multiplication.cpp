/* 
    matrix_multiplication.cpp
    butterfly matrix multiplication
*/
#include "matrix_multiplication.hpp"

void butterfly_matrix_multiplication (double *x,
                                      double *thetas,
                                      int half_len,
                                      int x_start_idx,
                                      int theta_start_idx,
                                      int mat_idx) {
    /*
       Compute butterfly matrix - vector multiplication
       Recursively compute based on mat_idx
       dimension should be the power of 2

       Input Parameters
       ----------------
       x :: input vector
       thetas :: array of thetas
       half_len :: effective half length of x (=x_end_idx-x_start_idx)
       x_start_idx :: start index of x
       theta_start_idx :: start index of theta
       mat_idx :: category of matrix index (mat_idx < log(half_len))
    */

    double old_x[2];
    double new_x[2];
    double theta, costheta, sintheta;

    int iter_idx;
    int x_idx;
    // mat_idx = 0 
    if (mat_idx == 0) {
        for (iter_idx=0; iter_idx<half_len; iter_idx++) {
            x_idx = x_start_idx + iter_idx;
            old_x[0] = x[x_idx];
            old_x[1] = x[x_idx+half_len];
            theta = thetas[theta_start_idx+iter_idx];
            costheta = cos(theta);
            sintheta = sin(theta);
            new_x[0] = costheta * old_x[0] - sintheta * old_x[1];
            new_x[1] = sintheta * old_x[0] + costheta * old_x[1];
            x[x_idx] = new_x[0];
            x[x_idx+half_len] = new_x[1];
        }
    }
    // mat_idx > 0
    else {
        butterfly_matrix_multiplication(x,
                                        thetas,
                                        half_len/2,
                                        x_start_idx,
                                        theta_start_idx,
                                        mat_idx-1);
        butterfly_matrix_multiplication(x,
                                        thetas,
                                        half_len/2,
                                        x_start_idx+half_len,
                                        theta_start_idx+(half_len/2),
                                        mat_idx-1);
    }
}
