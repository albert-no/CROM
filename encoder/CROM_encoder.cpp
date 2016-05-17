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

void CROM_encoder(double *x, int x_dim, int L, int *m_array, bool print_l2norm) {
    /*
       CROM encoder
       Assume k=1

       Input Parameters
       ______________-_
       x :: input vector
       x_dim :: dimension of x
       L :: number of iterations
       m_array :: array of massages
       print_l2norm :: whether printing intermediate l2 norm
    */
    double n = static_cast<double> (x_dim);
    double logn = log(n);
    int long_logn = static_cast<int> (logn);

    int half_mat_dim = x_dim/2;
    int x_start_idx = 0;
    int x_end_idx = x_dim;
    int theta_start_idx = 0;
    int theta_end_idx = half_mat_dim;
    int mat_idx;

    double *thetas = (double *)malloc (sizeof(double)*half_mat_dim);

    double scale = sqrt(n*(1-exp(-2*log(n)/n)));
    double scale_factor= exp(-log(n)/n);
    double uni_rand;

    // at i-th iterationof 
    // scale = [sqrt(n*(1-exp(-2*R/rawL))) * exp(-i*R/rawL)
    // R/L = log(n)/n
    int iter_idx;
    int x_cpy_idx;
    int theta_idx;
    int m;
    double l2norm;

    // Set random seed for thetas
    srand(THETA_SEED);
    fftw_plan p;
    p = fftw_plan_r2r_1d(x_dim, x, x, FFTW_REDFT10, FFTW_MEASURE);
    for (iter_idx=0; iter_idx<L; iter_idx++) {
        printf("iteration = %d\n", iter_idx);
        mat_idx = iter_idx % long_logn;

        // generate thetas from random seed
        for (theta_idx=0; theta_idx<half_mat_dim; theta_idx++) {
            uni_rand = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            thetas[theta_idx] = uni_rand * M_PI;
        }
        // multiply butterfly matrix
        butterfly_matrix_multiplication(x,
                                        thetas,
                                        half_mat_dim,
                                        x_start_idx,
                                        x_end_idx,
                                        theta_start_idx,
                                        theta_end_idx,
                                        mat_idx);
        // run dct2
        if (print_l2norm) {
            print_vector(x, x_dim);
        }
        fftw_execute(p);
        // normalize after dct2
        if (print_l2norm) {
            print_vector(x, x_dim);
        }
        normalize_vector(x, x_dim);
        scale *= scale_factor;

        // run CROM
        if (print_l2norm) {
            print_vector(x, x_dim);
        }
        m = CROM_step(x, x_dim, scale);
        m_array[iter_idx] = m;
        if (print_l2norm) {
            l2norm = compute_l2(x, x_dim);
            l2norm /= n;
            printf("m = %d, l2norm = %f\n", m, l2norm);
        }
    }
    fftw_destroy_plan(p);
    free(thetas);
}
