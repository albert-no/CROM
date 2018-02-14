// cromq_encoder.cpp
#include "cromq_encoder.hpp"
#include <iomanip>

using namespace Eigen;

CROMq_encoder::CROMq_encoder(std::string name_in,
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

    cov.resize(num_x, std::vector<double>(num_x));
    v_mat.resize(num_x, std::vector<double>(num_x));
    x_array.resize(x_dim, std::vector<double>(num_x));
}

CROMq_encoder::~CROMq_encoder() {
}

void CROMq_encoder::get_q_scores_and_mu(std::vector<std::vector<double>> &q_scores) {
    /* Read q_array from file
       convert char type to double */

    int row_idx = 0;
    int col_idx = 0;
    int subseq_idx;
    char qval;
    double value;

    std::ifstream file(fname);
    std::string line;

    if (verbose)
        std::cout << "Reading Q scores and computing mu" << std::endl;

    // initialize mu
    for (col_idx=0; col_idx<num_x; col_idx++)
        mu[col_idx] = 0;

    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        subseq_idx = 0;
        while(lineStream >> qval) {
            // convert charactor to double
            value = (double)qval;
            q_scores[row_idx][subseq_idx] = value;

            // update mu vector
            mu[subseq_idx] += value;
            subseq_idx++;
        }
        row_idx++;
    }
    // close file
    file.close();

    // normalize mu
    for (col_idx=0; col_idx<num_x; col_idx++) {
        mu[col_idx] /= (double) x_dim;
        if (verbose)
            std::cout << std::setprecision(3) << mu[col_idx] << " ";
    }
    if (verbose)
        std::cout << std::endl << std::endl;
}

void CROMq_encoder::get_cov(std::vector<std::vector<double>> &q_scores) {

    int row_idx, idx1, idx2;
    std::vector<double> q_temp(num_x);

    if (verbose) {
        std::cout << "Getting covariance matrix" << std::endl;
        std::cout << "num_x = " << num_x << std::endl;
    }
    // initialize cov
    for (idx1=0; idx1<num_x; idx1++)
        for(idx2=0; idx2<num_x; idx2++)
            cov[idx1][idx2] = 0;

    for (row_idx=0; row_idx<x_dim; row_idx++) {
        // subtract mu from q_scores first
        for (idx1=0; idx1<num_x; idx1++)
            q_temp[idx1] = q_scores[row_idx][idx1] - mu[idx1];

        // note: (estimated) covariance matrix is symmetric
        for (idx1=0; idx1<num_x; idx1++) {
            for(idx2=0; idx2<=idx1; idx2++) {
                cov[idx1][idx2] += (q_temp[idx1] * q_temp[idx2]);
            }
        }
    }

    // normalize cov
    for (idx1=0; idx1<num_x; idx1++) {
        cov[idx1][idx1] /= (double) x_dim;
        for(idx2=0; idx2<idx1; idx2++) {
            cov[idx1][idx2] /= (double) x_dim;
            cov[idx2][idx1] = cov[idx1][idx2];
        }
    }
    if (verbose) {
        for (idx1=0; idx1<num_x; idx1++) {
            for(idx2=0; idx2<num_x; idx2++) {
                std::cout << std::setprecision(3) << cov[idx1][idx2] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
}

void CROMq_encoder::do_svd() {
    int row_idx, col_idx;
    MatrixXf cov_eigen(num_x, num_x);
    MatrixXf sValues, vMatrix;

    if (verbose)
        std::cout << "Running SVD" << std::endl;

    for (row_idx=0; row_idx<num_x; row_idx++) {
        for (col_idx=0; col_idx<num_x; col_idx++) {
            cov_eigen(row_idx, col_idx) = cov[row_idx][col_idx];
        }
    }
    JacobiSVD<MatrixXf> svd(cov_eigen, ComputeThinU | ComputeThinV);

    sValues = svd.singularValues();
    vMatrix = svd.matrixV();
    for (row_idx=0; row_idx<num_x; row_idx++) {
        std_array[row_idx] = sqrt(sValues(row_idx));

        for (col_idx=0; col_idx<num_x; col_idx++) {
            v_mat[row_idx][col_idx] = vMatrix(row_idx, col_idx);
        }
    }
    if (verbose) {
        for (row_idx=0; row_idx<num_x; row_idx++) {
            std::cout << std_array[row_idx] << " ";
            for (col_idx=0; col_idx<num_x; col_idx++) {
                std::cout << std::setprecision(3) << v_mat[row_idx][col_idx] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }

    // write svd_params
    write_svd_params();
}

void CROMq_encoder::write_svd_params() {
    std::string mu_fname, v_mat_fname, std_array_fname;

    mu_fname = get_mu_fname(name, id);
    v_mat_fname = get_v_mat_fname(name, id);
    std_array_fname = get_std_array_fname(name, id);

    std::ofstream mu_file(mu_fname);
    std::ofstream v_mat_file(v_mat_fname);
    std::ofstream std_array_file(std_array_fname);

    int row_idx, col_idx;
    for (row_idx=0; row_idx<num_x; row_idx++) {
        mu_file << mu[row_idx] << std::endl;
        std_array_file << std_array[row_idx] << std::endl;
        for (col_idx=0; col_idx<num_x; col_idx++) {
            v_mat_file << v_mat[row_idx][col_idx] << " ";
        }
        v_mat_file << std::endl;
    }

    mu_file.close();
    std_array_file.close();
    v_mat_file.close();
}

void CROMq_encoder::normalize_q_scores(std::vector<std::vector<double>> &q_scores) {
    /* normalize q_scores with mu vector and V matrix*/
    int row_idx, col_idx, tmp_idx;
    if (verbose)
        std::cout << "Normalizing Q scroes" << std::endl;
    for (row_idx=0; row_idx<x_dim; row_idx++) {
        for(col_idx=0; col_idx<num_x; col_idx++) {
            x_array[row_idx][col_idx] = 0;
            // compute V^T * (Q-mu)
            for(tmp_idx=0; tmp_idx<num_x; tmp_idx++) {
                x_array[row_idx][col_idx] += (v_mat[tmp_idx][col_idx] * (q_scores[row_idx][tmp_idx] - mu[tmp_idx]));
            }
            // normalize with std_array
            x_array[row_idx][col_idx] /=  std_array[col_idx];
        }
    }
}

void CROMq_encoder::get_x_array() {
    /* Read q_scores from file
       Then do SVD and extract R_array, and std_array
       Then convert it to x_array
       x_array[i] will be an x_dim-dimensional array that contains i-th
       subsequence. Note that the x_array will be look like the transpose of
       input file. */

    std::vector<std::vector<double>> q_scores;
    int row_idx;
    int col_idx;

    //allocate memory
    q_scores.resize(x_dim, std::vector<double>(x_dim));

    // get q_scores and update mu vector
    get_q_scores_and_mu(q_scores);

    // get cov
    get_cov(q_scores);

    // do SVD
    do_svd();

    // normalized q score. this sets x_array
    normalize_q_scores(q_scores);

    if (verbose) {
        std::cout << std::endl;
        std::cout << "print x array" << std::endl;
        for (row_idx=0; row_idx<x_dim; row_idx++) {
            for (col_idx=0; col_idx<num_x; col_idx++) {
                std::cout << x_array[row_idx][col_idx] << " ";
            }
            std::cout << std::endl;
        }
    }
}

void CROMq_encoder::run() {

    std::string subname;
    int subseq_idx;

    // get_x_array
    get_x_array();

    //allocate_rate
    if (verbose)
        std::cout << "allocating rate" << std::endl;
    num_nonzero_rate = allocate_rate(std_array, R_array, num_x, rd_param, R_overall);

    // FIXME define default constructor
    // FIXME define copyable object and use vector
    // define fill member function in CROM_encoder
    // define nondefault constructor using fill member function
    // CROM_encoder* subenc;
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
