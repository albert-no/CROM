/*
   CROM_util.cpp
*/
#include "CROM_util.hpp"

int find_max_index(double *x, int x_dim) {
    /*
       Find the index of maximum element

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

double compute_l2(double *x, int x_dim) {
    /*
       Computes l2 norm of the vector

       Input Parameters
       ______________-_
       x :: input vector
       x_dim :: dimension of x

       Return Parameters
       _________________
       l2norm :: l2norm
    */

    int iter_idx;
    double l2norm = 0;

    for (iter_idx=0; iter_idx<x_dim; iter_idx++) {
        l2norm += (x[iter_idx] * x[iter_idx]);
    }
    return l2norm;
}

void print_vector(double *x, int x_dim) {
    /*
       Print vector

       Input Parameters
       ______________-_
       x :: input vector
       x_dim :: dimension of x

       Return Parameters
       _________________
    */

    int iter_idx;
    for (iter_idx=0; iter_idx<x_dim; iter_idx++) {
        printf("%6.3f, ", x[iter_idx]);
    }
    printf("\n\n");
}

void normalize_vector(double *x, int x_dim) {
    /*
       Normalize vector after dct2

       DCT II of fftw3 is a standard one. In order to normalize it
       (or make DCT-II matrix to be an orthogonal matrix), wwe need to scale
       X_0 by \sqrt{1/n} and all other by \sqrt{2/N}.

       Input Parameters
       ______________-_
       x :: input vector
       x_dim :: dimension of x

       Return Parameters
       _________________
    */

    int iter_idx;
    double n = static_cast <double> (x_dim);
    double inverse_sqrtn = sqrt(1/4.0/n);
    double inverse_halfsqrtn = sqrt(1/2.0/n);

    x[0] *= inverse_sqrtn;
    for (iter_idx=1; iter_idx<x_dim; iter_idx++) {
        x[iter_idx] *= inverse_halfsqrtn;
    }
}
