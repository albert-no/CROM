// CROMq_encoder.cpp
#include "CROMq_encoder.hpp"

CROMq_encoder::CROM_encoder(std::string name,
                            std::string fname_in,
                            int num_x_in,
                            int x_dim_in,
                            double *std_array_in,
                            double *R_array_in) {
    name = name_in;
    fname = fname_in;
    num_x = num_x_in;
    x_dim = x_dim_in;

    int sub_iter;
    R_array = new double[num_x];
    std_array = new double[num_x];
    for (sub_iter=0; sub_iter<num_x; sub_iter++) {
        std_array[sub_iter] = std_array_in[sub_iter];
        R_array[sub_iter] = R_array_in[sub_iter];
    }
}

CROMq_encoder::~CROMq_encoder() {
    if (R_array) {
        delete[] R_array;
    }
    if (std_array) {
        delete[] std_array;
    }
}

void CROMq_encoder::run() {
    /* read x_array from file
       x_array[i] will be an x_dim-dimensional array that contains i-th
       subsequence. Note that the x_array will be look like the transpose of
       input file. */
    double *x_array;
    x_array = new double[num_x][x_dim];
    std::ifstream file(fname);
    std::string line;
    int row_idx = 0;
    int subseq_idx;
    double value;
    while(std::getline(file, line)) {
        std::stringstream lineStream(line);
        subseq_idx = 0;
        while(lineStream >> value) {
            // normalize with std so that x_array contains N(0,1)
            x_array[subseq_idx][row_idx] = value / std_array[subseq_idx];
            subseq_idx++;
        }
        row_idx++;
    }

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

    // close file
    file.close();

    // clear the array
    delete[] x_array;
}
