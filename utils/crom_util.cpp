// crom_util.cpp
#include "crom_util.hpp"

int find_max_index(double *x, int x_dim) {
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
    int iter_idx;
    double l2norm = 0;

    for (iter_idx=0; iter_idx<x_dim; iter_idx++) {
        l2norm += (x[iter_idx] * x[iter_idx]);
    }
    return l2norm;
}

double compute_l2_dist(double *x, double *x_hat, int x_dim) {
    int iter_idx;
    double l2dist = 0;
    double diff;

    for (iter_idx=0; iter_idx<x_dim; iter_idx++) {
        diff = x_hat[iter_idx] - x[iter_idx];
        l2dist += (diff * diff);
    }
    return l2dist;
}

void print_vector(double *x, int x_dim) {
    int iter_idx;
    for (iter_idx=0; iter_idx<x_dim; iter_idx++) {
        printf("%6.3f, ", x[iter_idx]);
    }
    printf("\n\n");
}

void unnormalize_vector(double *x, int x_dim) {
    int iter_idx;
    double n = static_cast <double> (x_dim);
    double sqrtn = sqrt(1.0/n);
    double halfsqrtn = sqrt(1.0/2.0/n);

    x[0] = x[0] * sqrtn;
    for (iter_idx=1; iter_idx<x_dim; iter_idx++) {
        x[iter_idx] = x[iter_idx] * halfsqrtn;
    }
}

void copy_vector(double *x, double *xout, int x_dim) {
    int iter_idx;

    for (iter_idx=0; iter_idx<x_dim; iter_idx++) {
        x[iter_idx] = xout[iter_idx];
    }
}

void normalize_then_copy_vector(double *x, double *xout, int x_dim) {
    int iter_idx;
    double n = static_cast <double> (x_dim);
    double inverse_sqrtn = sqrt(1/4.0/n);
    double inverse_halfsqrtn = sqrt(1/2.0/n);

    x[0] = inverse_sqrtn * xout[0];
    for (iter_idx=1; iter_idx<x_dim; iter_idx++) {
        x[iter_idx] = inverse_halfsqrtn * xout[iter_idx];
    }
}

void generate_theta_from_seed(double *thetas, int theta_dim, int seed, bool sign) {
    int theta_idx;
    double uni_rand;
    double double_sign;

    double_sign = (sign) ? 1.0 : -1.0;

    srand(seed);
    for (theta_idx=0; theta_idx<theta_dim; theta_idx++) {
        uni_rand = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        thetas[theta_idx] = double_sign * uni_rand * M_PI;
    }
}