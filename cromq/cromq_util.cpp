// cromq_util.cpp
#include "cromq_util.hpp"


int allocate_rate(std::vector<double> &std_array, std::vector<double> &R_array, int num_x,
                  double rd_param, double R_overall) {
    /* Allocate rate according to the std_array
       Assuming D = e^{-1.4R} */

    std::vector<double> log_var_array(num_x);
    double sum_log_var = 0;
    double logD;
    int idx;
    int num_nonzero_rate;

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
    return num_nonzero_rate;
}

std::string get_mu_fname(std::string name, int id) {
    std::string mu_fname;
    mu_fname = "mu_id_" + std::to_string(id) + "_" + name + ".txt";
    return mu_fname;
}

std::string get_v_mat_fname(std::string name, int id) {
    std::string v_mat_fname;
    v_mat_fname = "v_mat_id_" + std::to_string(id) + "_" + name + ".txt";
    return v_mat_fname;
}

std::string get_std_array_fname(std::string name, int id) {
    std::string std_array_fname;
    std_array_fname = "std_array_id_" + std::to_string(id) + "_" + name + ".txt";
    return std_array_fname;
}

double compute_distortion(std::string ifname, std::string ofname, int num_x, int x_dim) {
    int row_idx, col_idx;
    char q_orig, q_recon;
    double distortion = 0;
    double diff;

    std::ifstream ifile(ifname);
    std::ifstream ofile(ofname);
    std::string iline, oline;

    for (row_idx=0; row_idx<x_dim; row_idx++) {
        std::getline(ifile, iline);
        std::getline(ofile, oline);
        std::stringstream ilineStream(iline);
        std::stringstream olineStream(oline);

        for (col_idx=0; col_idx<num_x; col_idx++) {
            ilineStream >> q_orig;
            olineStream >> q_recon;
            diff = ((double)q_orig - (double)q_recon);
            distortion += (diff*diff);
        }
    }
    distortion /= (num_x * x_dim);
    return distortion;
}
