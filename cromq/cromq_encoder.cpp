// cromq_encoder.cpp
#include "cromq_encoder.hpp"


CROMq_encoder::CROMq_encoder(std::string name_in,
                            std::string fname_in,
                            int id_in,
                            int num_x_in,
                            int x_dim_in,
                            double rd_param_in,
                            std::vector<double> R_array_in,
                            std::vector<double> mu_in,
                            std::vector<double> std_array_in,
                            std::vector<std::vector<double>> v_mat_in,
                            bool verbose_in) {
    int row_idx, col_idx;

    name = name_in;
    fname = fname_in;
    id = id_in;
    num_x = num_x_in;
    x_dim = x_dim_in;
    rd_param = rd_param_in;
    verbose = verbose_in;

    mu.resize(num_x);
    std_array.resize(num_x);
    v_mat.resize(num_x, std::vector<double>(num_x));
    R_array.resize(num_x);

    // copy svd params and allocated rates
    for (row_idx=0; row_idx<num_x; row_idx++) {
        mu[row_idx] = mu_in[row_idx];
        std_array[row_idx] = std_array_in[row_idx];
        R_array[row_idx] = R_array_in[row_idx];
        for (col_idx=0; col_idx<num_x; col_idx++) {
            v_mat[row_idx][col_idx] = v_mat_in[row_idx][col_idx];
        }
    }

    x_array.resize(x_dim, std::vector<double>(num_x));
}

CROMq_encoder::~CROMq_encoder() {
}

void CROMq_encoder::get_x_array() {
    /* Read q_scores from file
       Then do SVD and extract R_array, and std_array
       Then convert it to x_array
       x_array[i] will be an x_dim-dimensional array that contains i-th
       subsequence. Note that the x_array will be look like the transpose of
       input file.

       normalize q_scores with mu vector and V matrix */
    int row_idx = 0;
    int col_idx, tmp_idx;
    if (verbose)
        std::cout << "Get and Normalizing Q scores" << std::endl;

    std::ifstream file(fname);
    std::string line;

    while (std::getline(file, line)) {
        for (col_idx=0;col_idx<num_x;col_idx++) {
            x_array[row_idx][col_idx] = 0;
            // compute V^T * (Q-mu)
            for (tmp_idx=0; tmp_idx<num_x; tmp_idx++) {
                x_array[row_idx][col_idx] += (v_mat[tmp_idx][col_idx] * (line.at(tmp_idx) - mu[tmp_idx]));
            }
            // normalize with std_array
            x_array[row_idx][col_idx] /=  std_array[col_idx];
        }
        row_idx++;
    }
}

void CROMq_encoder::run() {

    std::string subname;
    int subseq_idx;

    // get_x_array
    get_x_array();

    // FIXME define copyable object and use vector
    // define fill member function in CROM_encoder
    // define nondefault constructor using fill member function
    // std::vector<CROM_encoder> subenc_array;

    for (subseq_idx=0; subseq_idx<num_x; subseq_idx++) {
        // declare sub_encoder
        subname = (name + "_id_" + std::to_string(id) + "_subid_" +
                   std::to_string(subseq_idx));
        if (verbose) {
            std::cout << "x_dim =  " << x_dim << std::endl;
        }
        CROM_encoder subenc(subname, x_dim, R_array[subseq_idx], false);

        // set x in sub_encoder
        if (verbose)
            std::cout << "Assigning normalized q scores to " << subseq_idx << "-th subencoder" << std::endl;
        subenc.set_x_from_array(x_array, subseq_idx, true);

        // run subencoder
        if (verbose) {
            std::cout << "Rate = " << R_array[subseq_idx] << std::endl;
            std::cout << "Running " << subseq_idx << "-th subencoder" << std::endl;
        }
        subenc.run();

        // write m_array in binary file
        subenc.write_m_array(true);
    }
}
