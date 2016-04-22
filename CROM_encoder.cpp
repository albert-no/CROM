/* 
   CROM_encoder.cpp
*/
#include "CROM_encoder.hpp"

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
    if (max_idx > BLOCKLENGTH) {
        cout << "ERROR!" << endl;
    }
    return max_idx;
}

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

void CROM_encoder(double *x, int x_dim, int L, int *m_array) {
    /*
       CROM encoder
       Assume k=1

       Input Parameters
       ______________-_
       x :: input vector
       x_dim :: dimension of x
       L :: number of iterations
       m_array :: array of massages
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

    double *thetas = (double *)malloc (half_mat_dim);
//    double *xtmp = (double *)malloc (BLOCKLENGTH);

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

    // Set random seed for thetas
    srand(THETA_SEED);
    fftw_plan p;
    p = fftw_plan_r2r_1d(BLOCKLENGTH, x, x, FFTW_REDFT10, FFTW_ESTIMATE);
    // p = fftw_plan_r2r_1d(BLOCKLENGTH, x, xtmp, FFTW_REDFT10, FFTW_ESTIMATE | FFTW_IN_PLACE);
    for (iter_idx=0; iter_idx<L; iter_idx++) {
        cout << "iteration = " << iter_idx << endl;
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
        fftw_execute(p);
        //for (x_cpy_idx=0; x_cpy_idx<BLOCKLENGTH; x_cpy_idx++) {
        //    x[x_cpy_idx] = xtmp[x_cpy_idx]/2.0;
        //}
        scale *= scale_factor;

        // run CROM
        m = CROM_step(x, x_dim, scale);
        m_array[iter_idx] = m;
    }
    fftw_destroy_plan(p);
    free(thetas);
    // free(xtmp);
}
