// cromq_encoder.cpp
#include "cromq_encoder.hpp"

using namespace Eigen;

CROMq_encoder::CROMq_encoder(std::string name_in,
                            std::string fname_in,
                            int num_x_in,
                            int x_dim_in,
                            double rd_param_in,
                            double R_overall_in) {
    int row_idx;

    name = name_in;
    fname = fname_in;
    num_x = num_x_in;
    x_dim = x_dim_in;
    rd_param = rd_param_in;
    R_overall= R_overall_in;

    std_array = new double[num_x];
    mu = new double[num_x];
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
    for (col_idx=0; col_idx<num_x; col_idx++)
        mu[col_idx] /= (double) x_dim;
}

void CROMq_encoder::get_cov(double** q_scores) {

    int row_idx, idx1, idx2;
    double* q_temp = new double[num_x];

    std::cout << "Getting covariance matrix" << std::endl;
    std::cout << "num_x = " << num_x << std::endl;
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
    for (idx1=0; idx1<num_x; idx1++) {
        cov[idx1][idx1] /= (double) x_dim;
        for(idx2=0; idx2<idx1; idx2++) {
            cov[idx1][idx2] /= (double) x_dim;
            cov[idx2][idx1] = cov[idx1][idx2];
        }
    }

    // free q_temp
    delete[] q_temp;
}

void CROMq_encoder::do_svd() {
    int row_idx, col_idx;
    MatrixXf cov_eigen(num_x, num_x);
    MatrixXf sValues, vMatrix;

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
        std_array[row_idx] = sValues(row_idx);
        for (col_idx=0; col_idx<num_x; col_idx++) {
            v_mat[row_idx][col_idx] = vMatrix(row_idx, col_idx);
        }
    }
}

void CROMq_encoder::normalize_q_scores(double** q_scores) {
    /* normalize q_scores with mu vector and V matrix*/
    int row_idx, col_idx, tmp_idx;
    std::cout << "Normalizing Q scroes" << std::endl;
    for (row_idx=0; row_idx<x_dim; row_idx++) {
        for(col_idx=0; col_idx<num_x; col_idx++) {
            x_array[row_idx][col_idx] = 0;
            // compute V * (Q-mu)
            for(tmp_idx=0; tmp_idx<num_x; tmp_idx++) {
                x_array[row_idx][col_idx] += v_mat[col_idx][tmp_idx] * (q_scores[row_idx][tmp_idx] - mu[tmp_idx]);
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

    // normalized q score. this sets x_array
    normalize_q_scores(q_scores);

    // free q_scores array
    for (row_idx=0; row_idx<x_dim; row_idx++)
        delete [] q_scores[row_idx];
    delete[] q_scores;
}

void CROMq_encoder::allocate_rate() {
    /* Allocate rate according to the std_array
       Assuming D = e^{-1.4R} */

    double* log_var_array = new double[num_x];
    double sum_log_var = 0;
    double logD;
    int idx;

    std::cout << "Allocating Rates" << std::endl;

    for (idx=0; idx<num_x; idx++) {
        log_var_array[idx] = 2 * log(std_array[idx]);
        sum_log_var += log_var_array[idx];
    }

    for (num_nonzero_rate=num_x; num_nonzero_rate>0; num_nonzero_rate--) {
        // solve :: (1/num_x) \sum log(sigma_i^2/D) = rd_param * R
        logD = sum_log_var / (double)num_x - rd_param * R_overall;
        if (logD <= log_var_array[num_nonzero_rate-1]) {
            for (idx=0; idx<num_nonzero_rate; idx++) {
                // allocate rates with following rule
                // R_i = (1/rd_param) * log(sigma^2/D)
                R_array[idx] = (log_var_array[idx] - logD) / rd_param;
            }
            break;
        }

        // assign zero rate for those subsequence with low variance
        R_array[num_nonzero_rate-1] = 0;
    }

    delete[] log_var_array;
}

void CROMq_encoder::run() {

    std::string subname;
    int subseq_idx;

    // get_x_array
    get_x_array();
    //allocate_rate
    allocate_rate();

    // FIXME define default constructor
    // define fill member function in CROM_encoder
    // define nondefault constructor using fill member function
    // CROM_encoder* subenc;
    std::vector<CROM_encoder> subenc_array;

    // subenc = new CROM_encoder[num_x];

    for (subseq_idx=0; subseq_idx<num_x; subseq_idx++) {
        // declare sub_encoder
        subname = (name + "_id_" + std::to_string(0) + "_subid_" +
                   std::to_string(subseq_idx));
        CROM_encoder subenc(subname, x_dim, R_array[subseq_idx], false);

        // set x in sub_encoder
        std::cout << "Assigning normalized q scores to " << subseq_idx << "-th subencoder" << std::endl;
        subenc.set_x_from_array(x_array, subseq_idx, true);

        //add sub_encoder to the array of sub_encoders
        subenc_array.push_back(subenc);

        // run subencoder
        std::cout << "Rate = " << R_array[subseq_idx] << std::endl;
        std::cout << "Running " << subseq_idx << "-th subencoder" << std::endl;
        subenc_array[subseq_idx].run();

        // write m_array in binary file
        subenc_array[subseq_idx].write_m_array(true);
    }
}
