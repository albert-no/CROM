/* 
   CROM_encoder.cpp
*/
#include "CROM_encoder.hpp"

// Constructor
CROM_encoder::CROM_encoder(std::string name_in, int x_dim_in, double R_in, bool verbose_in) {
    name = name_in;
    x_dim = x_dim_in;
    R = R_in;
    verbose = verbose_in;

    // compute the number of iterations
    double n = static_cast<double> (x_dim);
    L = static_cast<int> (n*R/log(n));

    // allocate memory
    m_array.resize(L);
    l2_array.resize(L);
    x.resize(x_dim);
}

// Destructor
CROM_encoder::~CROM_encoder() {
}

void CROM_encoder::set_x(std::vector<double> &x_in) {
    int x_iter;
    for (x_iter=0; x_iter<x_dim; x_iter++) {
        x[x_iter] = x_in[x_iter];
    }
}

void CROM_encoder::copy_x(std::vector<double> &x_copy) {
    int iter_idx;
    for (iter_idx=0; iter_idx<x_dim; iter_idx++) {
        x_copy[iter_idx] = x[iter_idx];
    }
}

int CROM_encoder::get_L() {
    return L;
}

void CROM_encoder::copy_m_array(std::vector<int> &m_array_copy) {
    int iter_idx;
    for (iter_idx=0; iter_idx<L; iter_idx++) {
        m_array_copy[iter_idx] = m_array[iter_idx];
    }
}

void CROM_encoder::copy_l2_array(std::vector<double> &l2_array_copy) {
    int iter_idx;
    for (iter_idx=0; iter_idx<L; iter_idx++) {
        l2_array_copy[iter_idx] = l2_array[iter_idx];
    }
}

int CROM_encoder::step(double scale) {
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
    double n = static_cast<double> (x_dim);
    double logn = log(n);
    int long_logn = static_cast<int> (logn);

    int half_len = x_dim/2;
    int x_start_idx = 0;
    int theta_start_idx = 0;
    int mat_idx;

    std::vector<double> thetas(half_len);
    std::vector<double> x_out(x_dim);

    double scale = sqrt(n*(1-exp(-2*log(n)/n)));
    double scale_factor= exp(-log(n)/n);
    double uni_rand;

    // At i-th iteration,
    // scale = sqrt(n*(1-exp(-2*R/rawL))) * exp(-i*R/rawL)
    // Note that R/L = log(n)/n
    int iter_idx;
    int theta_idx;
    int m;
    double l2norm;

    // Set random seed for thetas
    for (iter_idx=0; iter_idx<L; iter_idx++) {
        // before matrix multiplication
        if (verbose) {
            printf("iteration = %d\n", iter_idx);
            printf("Before matrix multiplication\n");
            print_vector(x, x_dim);
        }
        mat_idx = iter_idx % long_logn;

        // generate thetas for encoder (sign=true) with random seed=iter_idx
        generate_theta_from_seed(thetas, half_len, iter_idx, true);

        // multiply butterfly matrix
        butterfly_matrix_multiplication(x,
                                        thetas,
                                        half_len,
                                        x_start_idx,
                                        theta_start_idx,
                                        mat_idx);
        // run dct2
        FastDctFft::transform(x);

        // normalize after dct2
        normalize_vector(x, x_dim);

        // print after matrix multiplication
        if (verbose) {
            printf("After matrix multiplication\n");
            print_vector(x, x_dim);
        }
        // run CROM
        m = step(scale);
        m_array[iter_idx] = m;

        // print l2norm
        l2norm = compute_l2(x, x_dim);
        l2norm /= n;
        l2_array[iter_idx] = l2norm;
        if (verbose) {
            printf("Message and l2norm\n");
            printf("m = %d, l2norm = %f\n", m, l2norm);
        }

        // update scale with scale_factor
        scale *= scale_factor;
    }
}

void CROM_encoder::print_m_array() {
    int m_iter_idx;
    for (m_iter_idx=0; m_iter_idx<L; m_iter_idx++) {
        printf("%d\n", m_array[m_iter_idx]);
    }
}

void CROM_encoder::write_m_array(bool binary) {
    std::string filename;
    int line_idx;

    // if binary flag is on
    if (binary) {
        filename = "m_array_" + name + ".bin";
        std::ofstream m_outfile (filename, std::ios::binary);
        for (line_idx=0; line_idx<L; line_idx++) {
            m_outfile.write((char *)& m_array[line_idx], sizeof(m_array[line_idx]));
        }
        m_outfile.close();
    }
    else {
        filename = "m_array_" + name + ".txt";
        std::ofstream m_outfile;
        m_outfile.open(filename.c_str());
        for (line_idx=0; line_idx<L; line_idx++) {
            m_outfile << m_array[line_idx] << std::endl;
        }
        m_outfile.close();
    }
}
