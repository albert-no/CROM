// cromq_wrapper.cpp

#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

#include <Eigen/SVD>

#include "cromq_decoder.hpp"
#include "cromq_encoder.hpp"
#include "cromq_util.hpp"

using namespace Eigen;

class CROMq_wrapper
{

    std::string name;
    std::string svd_folder;
    std::string log_folder;

    int num_x;
    int xdim;

    long num_line;
    int num_file;

    double rd_param;
    double R_enc;

    bool verbose;

    int r_dec_num;
    std::vector<double> R_dec_array;
    std::vector<double> R_enc_per_col;

    std::vector<std::vector<double>> cov;

    std::vector<double> mu;
    std::vector<double> std_array;
    std::vector<std::vector<double>> v_mat;

    std::ofstream log_file;
    std::vector<std::string> subfnames;

    std::string mu_fname;
    std::string v_mat_fname;
    std::string std_array_fname;

    void split_qscores(bool skip);
    std::string get_subfname(int file_idx);
    void read_svd_params();
    void write_svd_params();
    void run_svd();
    void run_encoders();
    void run_decoders();
    void compute_distortions();

    public:
        CROMq_wrapper(std::string name_in,
                      int num_x_in,
                      int xdim_in,
                      double rd_param_in,
                      double R_enc_in,
                      int r_dec_num_in,
                      std::vector<double> R_dec_array_in,
                      std::string svd_folder_in,
                      std::string log_folder_in,
                      bool verbose_in);
        ~CROMq_wrapper();
        void run(bool encode, bool decode);
};
