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
#include "../cromq/cromq_util.hpp"


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


int main() {
    // Setup encoding parameters here
    // *********************
    std::string fname = "sample.qscore";
    int xdim = 4096;
    int num_x = 36;
    double rd_param = 1.4;
    double R_enc = 0.1;
    bool encode = true;
    bool decode = true;
    bool verbose = false;
    // *********************

    // Setup decoding parameters here
    // *********************
    int file_idx = 3;
    double R_dec = 0.05;
    // *********************
    
    // Generating log file
    std::string log_fname = "Renc_" + std::to_string(R_enc);
    log_fname += ("Rdec_" + std::to_string(R_dec));
    log_fname += ("_rd_" + std::to_string(rd_param));
    log_fname += ("_n_" + std::to_string(xdim));
    log_fname += ".log";

    std::ofstream log_file;
    log_file.open(log_fname);
    log_file << log_fname << std::endl;

    // check inputfile type
    if (fname.substr(fname.size()-7, fname.size()) != ".qscore") {
        std::cout << "Wrong input file name" << std::endl;
        exit(1);
    }
    std::string name = fname.substr(0, fname.size()-7);

    std::vector<std::string> subfnames;
    if (encode) {
        // split files
        file_idx = generate_subqscore_files(name, fname, subfnames, xdim);

        // XXX TBD take care of remainder (also when running cromq in the below)

        // Run CROMq encoder
        std::clock_t run_time;
        for (int cromq_idx=0; cromq_idx<file_idx-1; cromq_idx++) {
            run_time = std::clock();
            std::cout << "Processing " << subfnames[cromq_idx] << std::endl;
            CROMq_encoder enc(name,
                              subfnames[cromq_idx],
                              cromq_idx,
                              num_x,
                              xdim,
                              rd_param,
                              R_enc,
                              verbose);
            enc.run();
            run_time = std::clock() - run_time;
            double runtime_sec = ((double)run_time / (double)CLOCKS_PER_SEC);
            std::cout << "Enc took " << runtime_sec << " seconds" << std::endl << std::endl;
            log_file << "Enc took " << runtime_sec << " seconds" << std::endl << std::endl;
        }
    }
    if (decode) {
        // Run CROMq decoder
        std::clock_t run_time;
        if (!encode) { get_subfnames(name, subfnames, file_idx); }
        for (int cromq_idx=0; cromq_idx<file_idx-1; cromq_idx++) {
            run_time = std::clock();
            std::cout << "Processing " << subfnames[cromq_idx];
            std::cout << " at rate = " << R_dec << std::endl;
            CROMq_decoder dec(name,
                              subfnames[cromq_idx],
                              cromq_idx,
                              num_x,
                              xdim,
                              rd_param,
                              R_dec,
                              verbose);
            dec.run();
            run_time = std::clock() - run_time;
            double runtime_sec = ((double)run_time / (double)CLOCKS_PER_SEC);
            std::cout << "Dec took " << runtime_sec << " seconds" << std::endl << std::endl;
            log_file << "Dec took " << runtime_sec << " seconds" << std::endl << std::endl;
        }

        // Compute distortion
        double distortion;
        for (int cromq_idx; cromq_idx<file_idx-1; cromq_idx++) {
            std::string ofname = get_ofname(name, cromq_idx, R_dec);
            std::string ifname;
            std::stringstream ifname_tmp;
            ifname_tmp.str(std::string());
            ifname_tmp << std::setw(4) << std::setfill('0') << std::to_string(cromq_idx);
            ifname = "split_" + name + "_" + ifname_tmp.str() + ".subqscores";
            std::cout << ifname << " " << ofname << std::endl;
            distortion = compute_distortion(ifname, ofname, num_x, xdim);
            log_file << "Distortion (id " << cromq_idx << "): " << distortion << std::endl;
        }
    } 
    return 0;
}
