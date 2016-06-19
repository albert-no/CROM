#define TEST_BLOCKLENGTH 65536
#define TEST_HALFBLOCKLENGTH 32768

#include <cstdlib>
#include <fstream>

#include "../../encoder/CROM_encoder.hpp"
using namespace std;

int main() {
    double R = 0.1;
    double *x = new double[TEST_BLOCKLENGTH];
    int xdim = TEST_BLOCKLENGTH;
    bool verbose = false;
    ifstream x_infile;

    CROM_encoder enc (xdim, R, verbose);

    cout << "reading x input" << endl;
    string string_fname = "../input/x_input.txt";
    x_infile.open(string_fname.c_str());
    int read_line_idx;
    for (read_line_idx=0; read_line_idx<xdim; read_line_idx++) {
        x_infile >> x[read_line_idx];
    }
    x_infile.close();
    enc.set_x(x);
    delete[] x;

    cout << "Running CROM" << endl;
    enc.run();

    enc.print_m_array();
    return 0;
}
