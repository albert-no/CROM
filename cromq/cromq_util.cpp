// cromq_util.cpp
#include "cromq_util.hpp"

std::string get_log_fname(std::string logfolder, double R_enc, double rd_param, int xdim) {
    // Create log folder
    int status = mkdir(logfolder.c_str(), 0777);

    // Generating log file name
    std::stringstream logstream;
    logstream << logfolder << "/Renc" << std::fixed << std::setprecision(3) << R_enc;
    logstream << "_rd" << std::fixed << std::setprecision(3) << rd_param;
    logstream << "_n" << xdim;
    logstream << ".log";
    std::string log_fname = logstream.str();

    return log_fname;
}

int allocate_rate(std::vector<double> &std_array, std::vector<double> &R_array, int num_x,
                  double rd_param, double R_overall) {
    /* Allocate rate according to the std_array
       Assuming D = 2^{-1.4R} */

    std::vector<double> log_var_array(num_x);
    double sum_log_var = 0;
    double logD;
    int idx;
    int num_nonzero_rate;

    for (idx=0; idx<num_x; idx++) {
        log_var_array[idx] = 2 * log2(std_array[idx]);
        sum_log_var += log_var_array[idx];
    }

    for (num_nonzero_rate=num_x; num_nonzero_rate>0; num_nonzero_rate--) {
        // solve :: (1/num_x) \sum log2(sigma_i^2/D) = rd_param * R
        logD = (sum_log_var - num_x * rd_param * R_overall) / (double)num_nonzero_rate;
        if (logD <= log_var_array[num_nonzero_rate-1]) {
            for (idx=0; idx<num_nonzero_rate; idx++) {
                // allocate rates with following rule
                // R_i = (1/rd_param) * log2(sigma^2/D)
                R_array[idx] = (log_var_array[idx] - logD) / rd_param;
            }
            break;
        }

        // assign zero rate for those subsequence with low variance
        R_array[num_nonzero_rate-1] = 0;
        sum_log_var -= log_var_array[num_nonzero_rate-1];
    }
    return num_nonzero_rate;
}

std::string get_mu_fname(std::string name, int id) {
    std::string mu_fname;
    mu_fname = "svd_param/mu_id_" + std::to_string(id) + "_" + name + ".txt";
    return mu_fname;
}

std::string get_v_mat_fname(std::string name, int id) {
    std::string v_mat_fname;
    v_mat_fname = "svd_param/v_mat_id_" + std::to_string(id) + "_" + name + ".txt";
    return v_mat_fname;
}

std::string get_std_array_fname(std::string name, int id) {
    std::string std_array_fname;
    std_array_fname = "svd_param/std_array_id_" + std::to_string(id) + "_" + name + ".txt";
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

std::string get_ofname(std::string name, int id, double R_overall) {
    std::string ofname;
    ofname = (name + "_id_" + std::to_string(id) + "_Rate_" + std::to_string(R_overall) + "_out.txt");
    return ofname;
}

int generate_subqscore_files(std::string name,
                             std::string fname,
                             std::vector<std::string> &subfnames,
                             int xdim) {

    int line_idx = 0;
    int file_idx = 0;

    std::ifstream qscore_file(fname);
    std::string line;
    std::stringstream subfname_tmp;
    bool end_indicator = true;
    std::string subfname;
    while (end_indicator) {
        subfname_tmp.str(std::string());
        subfname_tmp << std::setw(4) << std::setfill('0') << file_idx;
        subfname = "split_" + name + "_" + subfname_tmp.str() + ".subqscores";
        subfnames.push_back(subfname);
        std::ofstream qscore_subfile(subfname);
        end_indicator = false;
        while(std::getline(qscore_file, line)) {
            qscore_subfile << line << std::endl;
            line_idx++;
            if (line_idx == xdim) {
                end_indicator = true;
                break;
            }
        }
        qscore_subfile.close();
        line_idx = 0;
        file_idx++;
    }
    qscore_file.close();

    return file_idx;
}


void get_subfnames(std::string name, std::vector<std::string> &subfnames, int file_idx) {
    std::stringstream subfname_tmp;
    std::string subfname;
    for (int iter_idx=0; iter_idx<file_idx; iter_idx++) {
        subfname_tmp.str(std::string());
        subfname_tmp << std::setw(4) << std::setfill('0') << iter_idx;
        subfname = "split_" + name + "_" + subfname_tmp.str() + ".subqscores";
        subfnames.push_back(subfname);
    }
}
