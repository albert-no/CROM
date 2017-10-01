// CROMq_encoder.cpp
#include "CROMq_encoder.hpp"

CROMq_encoder::CROMq_encoder(std::string name,
                            std::string fname_in,
                            int num_x_in,
                            int x_dim_in) {
    int row_idx;

    name = name_in;
    fname = fname_in;
    num_x = num_x_in;
    x_dim = x_dim_in;

    std_array = new double[num_x];
    mu_array = new double[num_x];
    R_array = new double[num_x];

    cov = new double*[num_x];
    for(row_idx=0; row_idx<num_x; row_idx++)
        cov[row_idx] = new double[num_x];

    v_mat = new double*[num_x];
    for(row_idx=0; row_idx<num_x; row_idx++)
        v_mat[row_idx] = new double[num_x];

    x_array = new double*[x_dim];
    for(row_idx=0; row_idx<x_dim; row_idx++)
        x_array[row_idx] = new double[num_x];
}

CROMq_encoder::~CROMq_encoder() {
    int row_idx;
    if (R_array)
        delete[] R_array;
    if (std_array)
        delete[] std_array;
    if (mu)
        delete[] mu;

    if (x_array) {
        for(row_idx=0; row_idx<x_dim; row_idx++)
            delete[] x_array[row_idx];
        delete[] x_array;
    }
    if (v_mat) {
        for (row_idx=0; row_idx<num_x; row_idx++)
            delete[] v_mat[row_idx];
        delete[] v_mat;
    }
    if (cov) {
        for (row_idx=0; row_idx<num_x; row_idx++)
            delete[] cov[row_idx];
        delete[] cov;
    }
}

void CROMq_encoder::allocate_rate() {
    /* Allocate rate according to the std_array

       Assuming D = e^{-1.4R} */
}

void CROMq_encoder::get_q_scores_and_mu(double** q_scores) {
    /* Read q_array from file
       convert char type to double */

    int row_idx = 0;
    int col_idx = 0;
    int subseq_idx;
    char qval;
    double value;

    std::ifstream file(fname);
    std::string line;

    // initialize mu
    for (col_idx=0; col_idx<num_x; col_idx++)
        mu[col_idx] = 0;

    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        subseq_idx = 0;
        while(lineStream >> qval) {
            // convert charactor to double
            value = (double)qval;
            q_array[row_idx][subseq_idx] = value;

            // update mu vector
            mu[subseq_idx] += value;
            subseq_idx++;
        }
        row_idx++;
    }
    // close file
    file.close();

    // normalize mu
    for (col_idx=0; col_idx<num_x; col_idx++)
        mu[col_idx] /= x_dim;
}

void CROMq_encoder::get_cov(double** q_scores) {

    int row_idx, idx1, idx2;
    double q_temp = new double[num_x];

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
                cov[idx1][idx2] += q_temp[idx1] * q_temp[idx2];
            }
        }
    }

    // normalize cov
    for (idx1=0; idx1<num_x; idx1++)
        for(idx2=0; idx2<num_x; idx2++)
            cov[idx1][idx2] /= x_dim;

    // free q_temp
    delete[] q_temp;
}

void CROMq_encoder::do_svd() {
    int row_idx, col_idx;
    MatrixXf cov_eigen;
    MatrixXf sValues, vMatrix;

    for (row_idx=0; row_idx<num_x; row_idx++) {
        for (col_idx=0; col_idx<num_x; col_idx++) {
            cov_eigen(row_idx, col_idx) = cov[row_idx][col_idx];
        }
    }
    JacobiSVD<MatrixXf> svd(cov_eigen, ComputeThinU | Compute ThinV);

    sValues = svd.singularValues();
    vMatrix = svd.matrixV();
    for (row_idx=0; row_idx<num_x; row_idx++) {
        std_array[row_idx] = sValues()(row_idx);
        for (col_idx=0; col_idx<num_x; col_idx++) {
            v_mat[row_idx][col_idx] = vMatrix(row_idx, col_idx);
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

    double** q_scores;
    int row_idx;

    //allocate memory
    q_scores = new double* [x_dim];
    for (row_idx=0; row_idx<x_dim; row_idx++)
        q_scores[row_idx] = new double[num_x];

    // get q_scores and update mu vector
    get_q_scores_and_mu(q_scores);

    // get cov
    get_cov(q_scores);

    // do SVD
    do_svd();

    // get x_array (normalized q score)
    get_x_array(q_scores);

    // free q_scores array
    for (row_idx=0; row_idx<x_dim; row_idx++)
        delete [] q_scores[row_idx];
    delete[] q_scores;
}


void CROMq_encoder::run() {

    std::string subname;
    for (subseq_idx=0; subseq_idx<num_x; subseq_idx++) {
        // declare sub_encoder
        // XXX TBD fix this
        subname = (name + "_id_" + std::tostring(0) + "_subid_" +
                   std::to_string(subseq_idx));
        CROM_encoder subenc(subname, x_dim, R_array[subseq_idx], false);
        // set x in sub_encoder
        subenc.set_x(x_array[subseq_idx]);
        // run subencoder
        subenc.run();
        // write m_array in binary file
        subenc.write_m_array(true);
    }
}
