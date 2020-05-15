// main.cpp

#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "../cromq/cromq_wrapper.hpp"
#include "../cromq/cromq_encoder.hpp"
#include "../cromq/cromq_util.hpp"


int main() {
    // Setup wrapper parameters here
    // *********************
    std::string name = "small";
    int num_x = 36;
    int xdim = 1024;
    double rd_param = 1.4;
    double R_enc = 0.2;
    int r_dec_num = 3;
    std::vector<double> R_dec_array = {0.05, 0.1, 0.2};
    bool verbose = false;
    std::string svd_folder = "svd_params";
    std::string log_folder = "logs";

    CROMq_wrapper cromq_wrapper(name,
                                num_x,
                                xdim,
                                rd_param,
                                R_enc,
                                r_dec_num,
                                R_dec_array,
                                svd_folder,
                                log_folder,
                                verbose);

    bool encode = true;
    bool decode = true;
    cromq_wrapper.run(encode, decode); 
    return 0;
}
