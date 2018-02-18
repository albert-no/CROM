// main.cpp
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "../cromq/cromq_decoder.hpp"
#include "../cromq/cromq_encoder.hpp"

int main() {
    // Setup parameters here
    // *********************
    std::string fname = "SRR494099.qscore";
    int xdim = 65536;
    int num_x = 36;
    double rd_param = 1.4;
    bool verbose = false;
    // double R_enc = 2;
    double R_enc = 2;
    int num_dec_pts = 4;
    // std::vector<double> R_dec = {0.5, 1, 1.5, 2};
    // *********************

    if (fname.substr(fname.size()-7, fname.size()) != ".qscore") {
        std::cout << "Wrong input file name" << std::endl;
        exit(1);
    }

    std::string name = fname.substr(0, fname.size()-7);

    // split files
    std::ifstream qscore_file(fname);
    std::string line;
    std::stringstream subfname_tmp;

    std::vector<std::string> subfnames;

    int line_idx = 0;
    int file_idx = 0;
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

    // XXX TBD take care of remainder (also when running cromq in the below)

    // Run CROMq
    std::clock_t run_time;
    for (int cromq_idx=0; cromq_idx<file_idx-1; cromq_idx++) {
        run_time = std::clock();
        std::cout << "Processing " << subfnames[cromq_idx] << std::endl;
        CROMq_encoder enc(name, subfnames[cromq_idx], cromq_idx, num_x, xdim, rd_param, R_enc, verbose);
        enc.run();
        run_time = std::clock() - run_time;
        std::cout << "It took " << ((double)run_time / (double)CLOCKS_PER_SEC) << " seconds" << std::endl << std::endl;
    }

    return 0;
}
