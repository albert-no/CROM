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


int generate_subqscore_files(std::string name, std::string fname, std::vector<std::string> &subfnames, int xdim) {
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

    return file_idx-1;
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


int main() {
    // Setup encoding parameters here
    // *********************
    std::string fname = "SRR494099.qscore";
    int xdim = 65536;
    int num_x = 36;
    double rd_param = 1.4;
    double R_enc = 0.1;
    char mode = 'd';
    bool verbose = false;
    // *********************

    // Setup decoding parameters here
    // *********************
    int file_idx = 2;
    int num_dec_pts = 2;
    std::vector<double> R_dec = {0.05, 0.09};
    // *********************

    if (fname.substr(fname.size()-7, fname.size()) != ".qscore") {
        std::cout << "Wrong input file name" << std::endl;
        exit(1);
    }
    std::string name = fname.substr(0, fname.size()-7);

    switch(mode) {
        case 'e': {
            // split files
            std::vector<std::string> subfnames;
            file_idx = generate_subqscore_files(name, fname, subfnames, xdim);

            // XXX TBD take care of remainder (also when running cromq in the below)

            // Run CROMq encoder
            std::clock_t run_time;
            for (int cromq_idx=0; cromq_idx<file_idx; cromq_idx++) {
                run_time = std::clock();
                std::cout << "Processing " << subfnames[cromq_idx] << std::endl;
                CROMq_encoder enc(name, subfnames[cromq_idx], cromq_idx, num_x, xdim, rd_param, R_enc, verbose);
                enc.run();
                run_time = std::clock() - run_time;
                std::cout << "It took " << ((double)run_time / (double)CLOCKS_PER_SEC) << " seconds" << std::endl << std::endl;
            }
        } break;

        case 'd': {
            // Run CROMq decoder
            std::clock_t run_time;
            std::vector<std::string> subfnames;
            get_subfnames(name, subfnames, file_idx);
            for (int cromq_idx=0; cromq_idx<file_idx; cromq_idx++) {
                for (int dec_idx=0; dec_idx<num_dec_pts; dec_idx++) {
                    run_time = std::clock();
                    std::cout << "Processing " << subfnames[cromq_idx] << " at rate = " << R_dec[dec_idx] << std::endl;
                    CROMq_decoder dec(name, subfnames[cromq_idx], cromq_idx, num_x, xdim, rd_param, R_dec[dec_idx], verbose);
                    dec.run();
                    run_time = std::clock() - run_time;
                    std::cout << "It took " << ((double)run_time / (double)CLOCKS_PER_SEC) << " seconds" << std::endl << std::endl;
                }
            }
        } break;

        default: {
            std::cout << "Wrong mode (it should be either `e` (encoding) or `d` (decoding)" << std::endl;
            exit(1);
        }
    }

    return 0;
}
