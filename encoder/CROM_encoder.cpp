/* 
   CROM_encoder.cpp
*/
#include "CROM_encoder.hpp"

int CROM_step(double *x, int x_dim, double scale) {
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
    max_idx = find_max_index(x, x_dim);
    for (iter_idx=0; iter_idx<x_dim; iter_idx++) {
        x[iter_idx] -= offset;
    }
    x[max_idx] += offset;
    x[max_idx] -= (sqrt((n-1.0)/n)*scale);
    return max_idx;
}

void CROM_encoder(double *x, int x_dim, int L, int *m_array, bool verbose) {
    /*
       CROM encoder
       Assume k=1

       Input Parameters
       ______________-_
       x :: input vector
       x_dim :: dimension of x
       L :: number of iterations
       m_array :: array of massages
       verbose :: whether printing intermediate l2 norm
    */
    double n = static_cast<double> (x_dim);
    double logn = log(n);
    int long_logn = static_cast<int> (logn);

    int half_len = x_dim/2;
    int x_start_idx = 0;
    int theta_start_idx = 0;
    int mat_idx;

    double *thetas = (double *)malloc (sizeof(double)*half_len);
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

    // Set random seed for thetas
    fftw_plan p;
    p = fftw_plan_r2r_1d(x_dim, x, x_out, FFTW_REDFT10, FFTW_MEASURE);
    for (iter_idx=0; iter_idx<L; iter_idx++) {
        printf("iteration = %d\n", iter_idx);
        mat_idx = iter_idx % long_logn;

        // generate thetas from random seed
        srand(iter_idx);
        for (theta_idx=0; theta_idx<half_len; theta_idx++) {
            uni_rand = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            thetas[theta_idx] = uni_rand * M_PI;
        }

        // multiply butterfly matrix
        butterfly_matrix_multiplication(x,
                                        thetas,
                                        half_len,
                                        x_start_idx,
                                        theta_start_idx,
                                        mat_idx);
        // run dct2
        fftw_execute(p);

        // normalize after dct2
        normalize_then_copy_vector(x, x_out, x_dim);

        // run CROM
        if (verbose) {
            print_vector(x, x_dim);
        }
        m = CROM_step(x, x_dim, scale);
        m_array[iter_idx] = m;

        // print l2norm
        l2norm = compute_l2(x, x_dim);
        l2norm /= n;
        printf("m = %d, l2norm = %f\n", m, l2norm);
        // update scale with scale_factor
        scale *= scale_factor;
    }
    fftw_destroy_plan(p);
    free(thetas);
    free(x_out);
}
