/* 
   CROM_encoder.cpp
*/
#include "CROM_encoder.hpp"

CROM_encoder::CROM_encoder(int x_dim_in, double R_in, bool verbose_in) {
    // Constructor
    x_dim = x_dim_in;
    R = R_in;
    verbose = verbose_in;

    // compute the number of iterations
    double n = static_cast<double> (x_dim);
    L = static_cast<int> (n*R/log(n));

    // allocate memory
    m_array = new int[L+1];
    x = new double[x_dim];
}


CROM_encoder::~CROM_encoder() {
    if (x)
        delete[] x;
    if (m_array)
        delete[] m_array;
}


void CROM_encoder::read_x(std::string filename) {
    std::ifstream x_infile;
    x_infile.open(filename.c_str());
    int read_line_idx;
    for (read_line_idx=0; read_line_idx<x_dim; read_line_idx++) {
        x_infile >> x[read_line_idx];
    }
    x_infile.close();
}


int CROM_encoder::step(double scale) {
    /*
       Single iteration of CROM encoder
       Assume k=1

       Parameters
       ----------
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


void CROM_encoder::run() {
    /* Run CROM encoder with k=1

    */
    double n = static_cast<double> (x_dim);
    double logn = log(n);
    int long_logn = static_cast<int> (logn);

    int half_len = x_dim/2;
    int x_start_idx = 0;
    int theta_start_idx = 0;
    int mat_idx;

    double *thetas = new double[half_len];
    double *x_out = new double[x_dim];

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
        m = step(scale);
        m_array[iter_idx] = m;

        // print l2norm
        l2norm = compute_l2(x, x_dim);
        l2norm /= n;
        printf("m = %d, l2norm = %f\n", m, l2norm);
        // update scale with scale_factor
        scale *= scale_factor;
    }
    fftw_destroy_plan(p);
    delete[] thetas;
    delete[] x_out;
}


void CROM_encoder::print_m_array() {
    int m_iter_idx;
    for (m_iter_idx=0; m_iter_idx<L; m_iter_idx++) {
        printf("%d\n", m_array[m_iter_idx]);
    }
}
