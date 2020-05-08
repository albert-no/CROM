/* 
   crom_decoder.cpp
*/
#include "crom_decoder.hpp"

// Constructor
CROM_decoder::CROM_decoder(std::string name_in, int x_dim_in, int L_in, bool verbose_in) {
    int x_iter;
    name = name_in;
    x_dim = x_dim_in;
    L = L_in;
    verbose = verbose_in;

    x_hat.resize(x_dim);
    for (x_iter=0; x_iter<x_dim; x_iter++) {
        x_hat[x_iter] = 0;
    }

    m_array.resize(L);
}

// Destructor
CROM_decoder::~CROM_decoder() {
}

void CROM_decoder::set_x_to_array(std::vector<std::vector<double>> &x_array, int idx, bool vertical) {
    int x_iter;
    if (vertical) {
        for (x_iter=0; x_iter<x_dim; x_iter++) {
            x_array[x_iter][idx] = x_hat[x_iter];
        }
    }
    else {
        for (x_iter=0; x_iter<x_dim; x_iter++) {
            x_array[idx][x_iter] = x_hat[x_iter];
        }
    }
}

void CROM_decoder::set_m_array(std::vector<int> &m_array_in) {
    int m_iter;
    for (m_iter=0; m_iter<L; m_iter++) {
        m_array[m_iter] = m_array_in[m_iter];
    }
}

void CROM_decoder::read_m_array(bool binary) {
    std::string filename;
    int line_idx;

    // if binary flag is on
    if (binary) {
        filename = "m_array_" + name + ".bin";
        std::ifstream m_infile (filename, std::ios::binary);

        int logn = int(log2(x_dim));
        int m_byte = 0;
        int m = 0;
        int line_idx = 0;
        int current_bit = 0;
        while (1) {
            m_byte = m_infile.get();
            if (m_byte == EOF) break;
            // 8 is for byte
            for (int m_idx=0; m_idx<8; m_idx++) {
                if (m_byte % 2 == 1) {
                    m |= (1<<current_bit);
                }
                current_bit++;
                m_byte = m_byte >> 1;
                if (current_bit == logn) {
                    m_array[line_idx] = m;
                    current_bit = 0;
                    m = 0;
                    line_idx ++;
                    if(line_idx == L) {
                        m_infile.close();
                        return;
                    }
                }
            }
            
        }
        m_infile.close();
    }
    else {
        filename = "m_array_" + name + ".txt";
        std::ifstream m_infile (filename);
        for (line_idx=0; line_idx<L; line_idx++) {
            m_infile >> m_array[line_idx];
        }
        m_infile.close();
    }
}

void CROM_decoder::copy_x_hat(std::vector<double> &x_hat_copy) {
    int x_iter;
    for (x_iter=0; x_iter<x_dim; x_iter++) {
        x_hat_copy[x_iter] = x_hat[x_iter];
    }
}

void CROM_decoder::step(double scale, int m) {
    int max_idx;
    int iter_idx;
    double n = static_cast<double> (x_dim);
    double offset;

    offset = -sqrt(1.0/n/(n-1.0)) * scale;
    for (iter_idx=0; iter_idx<x_dim; iter_idx++) {
        x_hat[iter_idx] += offset;
    }
    x_hat[m] -= offset;
    x_hat[m] += (sqrt((n-1.0)/n)*scale);
}

void CROM_decoder::run() {
    double n = static_cast<double> (x_dim);
    double log2n = log2(n);
    int long_logn = static_cast<int> (log2n);

    const int half_len = x_dim/2;
    int x_start_idx = 0;
    int theta_start_idx = 0;
    int mat_idx;
    int x_hat_idx;

    const int xout_size = x_dim;
    std::vector<double> xout(xout_size);

    std::vector<double> thetas_inv(half_len);

    double scale = sqrt(n*(1-exp(-2*log2n/n)));
    double scale_factor= exp(-log2n/n);
    double uni_rand;

    // at i-th iterationof 
    // scale = [sqrt(n*(1-exp(-2*R/rawL))) * exp(-i*R/rawL)
    // R/L = log2(n)/n
    int iter_idx;
    int theta_idx;
    int m;
    double l2norm;

    // setup scale
    scale = sqrt(n*(1-exp(-2*log2n/n))) * exp(-(L-1)*log2n/n);

    // Set random seed for thetas
    fftw_plan p;
    p = fftw_plan_r2r_1d(x_dim, x_hat.data(), xout.data(), FFTW_REDFT01, FFTW_MEASURE);
    for (iter_idx=L-1; iter_idx>=0; iter_idx--) {
        mat_idx = iter_idx % long_logn;

        // Decoding step
        m = m_array[iter_idx];
        step(scale, m);

        // unnormalize before idct2
        unnormalize_vector(x_hat, x_dim);

        // run idct2
        fftw_execute(p);

        // copy x_hat from xout
        copy_vector(x_hat, xout, x_dim);

        // generate thetas for decoder (sign=false) from random seed=iter_idx
        generate_theta_from_seed(thetas_inv, half_len, iter_idx, false);

        // multiply inverse butterfly matrix
        butterfly_matrix_multiplication(x_hat,
                                        thetas_inv,
                                        half_len,
                                        x_start_idx,
                                        theta_start_idx,
                                        mat_idx);
        // update scale with scale factor
        scale /= scale_factor;
    }
    fftw_destroy_plan(p);
}

void CROM_decoder::write_x_hat() {
    std::ofstream x_hat_outfile;
    std::string filename;

    filename = "x_hat_"+ name + ".txt";
    x_hat_outfile.open(filename.c_str());
    int line_idx;
    for (line_idx=0; line_idx<x_dim; line_idx++) {
        x_hat_outfile << x_hat[line_idx] << std::endl;
    }
    x_hat_outfile.close();
}
