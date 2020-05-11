// main.cpp

#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "../cromq/cromq_decoder.hpp"
#include "../cromq/cromq_encoder.hpp"
#include "../cromq/cromq_util.hpp"


int main() {
    // Setup encoding parameters here
    // *********************
    std::string fname = "SRR494099.qscore";
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
    double R_dec = 0.1;
    // *********************
    
    // Generating log file
    int status = mkdir("logs", 0777);
    std::stringstream logstream;
    logstream << "logs/Renc" << std::fixed << std::setprecision(3) << R_enc;
    logstream << "_Rdec" << std::fixed << std::setprecision(3) << R_dec;
    logstream << "_rd" << std::fixed << std::setprecision(3) << rd_param;
    logstream << "_n" << xdim;
    logstream << ".log";
    std::string log_fname = logstream.str();

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
        double enc_runtime = 0;
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
            std::cout << "Enc took " << runtime_sec << " seconds" << std::endl;
            log_file << "Enc took " << runtime_sec << " seconds" << std::endl;
            enc_runtime += runtime_sec;
        }
        std::cout << "Total: Enc took " << enc_runtime << " seconds" << std::endl << std::endl;
        log_file << "Total: Enc took " << enc_runtime << " seconds" << std::endl << std::endl;
    }
    if (decode) {
        double dec_runtime = 0;
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
            dec_runtime += runtime_sec;
        }
        std::cout << "Total: Dec took " << dec_runtime << " seconds" << std::endl;
        log_file << "Total: Dec took " << dec_runtime << " seconds" << std::endl;
    } 

    // Compute distortion
    std::cout << "Compute Distortion" << std::endl;
    double distortion_avg = 0;
    for (int cromq_idx=0; cromq_idx<file_idx-1; cromq_idx++) {
        double distortion;
        std::string ofname = get_ofname(name, cromq_idx, R_dec);
        std::string ifname;
        std::stringstream ifname_tmp;
        ifname_tmp << "split_" << name << "_";
        ifname_tmp << std::setfill('0') << std::setw(4) << cromq_idx;
        ifname_tmp << ".subqscores";
        ifname = ifname_tmp.str();
        distortion = compute_distortion(ifname, ofname, num_x, xdim);
        std::cout << "Distortion (id " << cromq_idx << "): " << distortion << std::endl;
        log_file << "Distortion (id " << cromq_idx << "): " << distortion << std::endl;
        distortion_avg += distortion;
    }
    distortion_avg /= (double)(file_idx-1);
    std::cout << "Average distortion : " << distortion_avg << std::endl;
    log_file << "Average distortion : " << distortion_avg << std::endl;
    return 0;
}
