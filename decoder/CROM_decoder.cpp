/* 
   CROM_decoder.cpp
*/
#include "CROM_decoder.hpp"

void CROM_decoding_step(double *xhat, int x_dim, double scale, int m) {
    /*
       Single iteration of CROM encoder
       Assume k=1

       Input Parameters
       ______________-_
       x :: input vector
       x_dim :: dimension of x
       scale :: scale factor of iteration
    */
    int max_idx;
    int iter_idx;
    double n = static_cast<double> (x_dim);
    double offset;

    offset = -sqrt(1.0/n/(n-1.0)) * scale;
    for (iter_idx=0; iter_idx<x_dim; iter_idx++) {
        xhat[iter_idx] += offset;
    }
    xhat[m] -= offset;
    xhat[m] += (sqrt((n-1.0)/n)*scale);
}

void CROM_decoder(double *xhat, int x_dim, int L, int *m_array, bool verbose) {
    /*
       CROM encoder
       Assume k=1

       Input Parameters
       ______________-_
       xhat :: reconstruction vector
               initialized with 0
       x_dim :: dimension of x
       L :: number of iterations
       m_array :: array of massages
    */
    double n = static_cast<double> (x_dim);
    double logn = log(n);
    int long_logn = static_cast<int> (logn);

    int half_len = x_dim/2;
    int x_start_idx = 0;
    int theta_start_idx = 0;
    int mat_idx;

    double *thetas_inv = (double *)malloc (sizeof(double)*half_len);
    double *x_out = (double *)malloc (sizeof(double)*x_dim);

    double scale = sqrt(n*(1-exp(-2*log(n)/n)));
    double scale_factor= exp(-log(n)/n);
    double uni_rand;

    // at i-th iterationof 
    // scale = [sqrt(n*(1-exp(-2*R/rawL))) * exp(-i*R/rawL)
    // R/L = log(n)/n
    int iter_idx;
    int theta_idx;
    int m;
    double l2norm;

    // setup scale
    scale = sqrt(n*(1-exp(-2*log(n)/n))) * exp(-(L)*log(n)/n);

    // Set random seed for thetas
    fftw_plan p;
    p = fftw_plan_r2r_1d(x_dim, xhat, x_out, FFTW_REDFT01, FFTW_MEASURE);
    for (iter_idx=L-1; iter_idx>=0; iter_idx--) {
        printf("iteration = %d\n", iter_idx);
        mat_idx = iter_idx % long_logn;

        // Decoding step
        m = m_array[iter_idx];
        CROM_decoding_step(xhat, x_dim, scale, m);

        // unnormalize before idct2
        unnormalize_vector(xhat, x_dim);
        // run idct2
        fftw_execute(p);
        // copy x from xout
        copy_vector(xhat, x_out, x_dim);

        scale /= scale_factor;

        // generate thetas from random seed
        srand(iter_idx);
        for (theta_idx=0; theta_idx<half_len; theta_idx++) {
            uni_rand = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            thetas_inv[theta_idx] = -uni_rand * M_PI;
        }
        // multiply inverse butterfly matrix
        butterfly_matrix_multiplication(xhat,
                                        thetas_inv,
                                        half_len,
                                        x_start_idx,
                                        theta_start_idx,
                                        mat_idx);
        if (verbose) {
            print_vector(xhat, x_dim);
        }
    }
    fftw_destroy_plan(p);
    free(thetas_inv);
    free(x_out);
}
