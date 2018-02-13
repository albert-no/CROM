// cromq_decoder.cpp
#include "cromq_decoder.hpp"
#include <iomanip>

using namespace Eigen;

CROMq_decoder::CROMq_decoder(std::string name_in,
                            std::string fname_in,
                            int id_in,
                            int num_x_in,
                            int x_dim_in,
                            double rd_param_in,
                            double R_overall_in,
                            bool verbose_in) {
    int row_idx;

    name = name_in;
    fname = fname_in;
    id = id_in;
    num_x = num_x_in;
    x_dim = x_dim_in;
    rd_param = rd_param_in;
    R_overall = R_overall_in;
    verbose = verbose_in;


    mu.resize(num_x);
    std_array.resize(num_x);
    R_array.resize(num_x);

    v_mat.resize(num_x, std::vector<double>(num_x));
    x_hat_array.resize(x_dim, std::vector<double>(num_x));
    q_scores.resize(x_dim, std::vector<char>(num_x));
}

CROMq_decoder::~CROMq_decoder() {
}

void CROMq_decoder::get_svd_params() {
    /* Read mu, v_mat, std_array from file
       convert char type to double */

    int row_idx;
    int col_idx;

    std::string mu_fname, v_mat_fname, std_array_fname;

    mu_fname = get_mu_fname(name, id);
    std_array_fname = get_std_array_fname(name, id);
    v_mat_fname = get_v_mat_fname(name, id);

    std::ifstream mu_file(mu_fname);
    std::ifstream std_array_file(std_array_fname);
    std::ifstream v_mat_file(v_mat_fname);

    if (verbose)
        std::cout << "Reading mu from " << mu_fname << std::endl;

    for (row_idx=0; row_idx<num_x; row_idx++) {
        mu_file >> mu[row_idx];
    }
    // close file
    mu_file.close();

    if (verbose)
        std::cout << "Reading std_array from " << std_array_fname << std::endl;

    for (row_idx=0; row_idx<num_x; row_idx++) {
        std_array_file >> std_array[row_idx];
    }
    // close file
    std_array_file.close();

    if (verbose)
        std::cout << "Reading v_mat from " << v_mat_fname << std::endl;

    row_idx = 0;
    for (row_idx=0; row_idx<num_x; row_idx++) {
        for (col_idx=0; col_idx<num_x; col_idx++) {
            v_mat_file >> v_mat[row_idx][col_idx];
        }
    }
    // close file
    v_mat_file.close();
}

char CROMq_decoder::q_score_round(double q_score) {
    int rounded;
    if (q_score < 33) {
        rounded = 33;
    }
    else if (q_score > 73) {
        rounded = 73;
    }
    else {
        rounded = round(q_score);
    }
    return rounded;
}

void CROMq_decoder::unnormalize_q_scores() {
    /* unnormalize q_scores with mu vector and V matrix
       Note that when we normalize, we compute X = V^T * (Q-mu), and then normalize with std_array */
    int row_idx, col_idx, tmp_idx;
    double q_score_temp;
    if (verbose)
        std::cout << "Unnormalizing Q scroes" << std::endl;

    for (row_idx=0; row_idx<x_dim; row_idx++) {
        for(col_idx=0; col_idx<num_x; col_idx++) {
            // unnormalize with std_array
            x_hat_array[row_idx][col_idx] *=  std_array[col_idx];
            // compute Q = V * X + mu
            q_score_temp = 0;
            q_scores[row_idx][col_idx] = 0;
            for(tmp_idx=0; tmp_idx<num_x; tmp_idx++) {
                q_score_temp += (v_mat[col_idx][tmp_idx] * x_hat_array[row_idx][tmp_idx]);
            }
            q_score_temp += mu[col_idx];
            q_scores[row_idx][col_idx] = q_score_round(q_score_temp);
        }
    }
}

void CROMq_decoder::write_q_scores() {
    std::string q_scores_fname = get_ofname(name, id, R_overall);
    std::ofstream q_scores_file(q_scores_fname);

    int row_idx, col_idx;
    for (row_idx=0; row_idx<x_dim; row_idx++) {
        for (col_idx=0; col_idx<num_x; col_idx++) {
            q_scores_file << q_scores[row_idx][col_idx];
        }
        q_scores_file << std::endl;
    }
    q_scores_file.close();
}

void CROMq_decoder::run() {

    std::string subname;

    int subseq_idx;
    int L;

    double n = static_cast<double> (x_dim);

    // set false if m-array file is non-binary
    bool binary = true;

    get_svd_params();

    //allocate_rate
    num_nonzero_rate = allocate_rate(std_array, R_array, num_x, rd_param, R_overall);

    if (verbose) {
        std::cout << "num_nonzero_rate = " << num_nonzero_rate << std::endl;
        for (subseq_idx=0; subseq_idx<num_x; subseq_idx++) 
            std::cout << R_array[subseq_idx] << " ";
        std::cout << std::endl;
    }

    for (subseq_idx=0; subseq_idx<num_x; subseq_idx++) {
        // declare sub_encoder
        subname = (name + "_id_" + std::to_string(id) + "_subid_" +
                   std::to_string(subseq_idx));

        L = static_cast<int> (n * R_array[subseq_idx] / log(n));
        if (verbose)
            std::cout << "L = " << L << std::endl;
        CROM_decoder subdec(subname, x_dim, L, false);

        // read m_array from file
        subdec.read_m_array(binary);

        // run subdecoder
        if (verbose) {
            std::cout << "Rate = " << R_array[subseq_idx] << std::endl;
            std::cout << "Running " << subseq_idx << "-th subdecoder" << std::endl;
        }
        subdec.run();

        // assign x_hat to x_hat_array 
        subdec.set_x_to_array(x_hat_array, subseq_idx, true);
    }
    // unnormalize q_scores
    unnormalize_q_scores();

    // write q_hat
    write_q_scores();
}
