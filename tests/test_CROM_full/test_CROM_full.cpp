#include <fstream>
#define TEST_BLOCKLENGTH 65536
#define TEST_HALFBLOCKLENGTH 32768

#include "../../encoder/CROM_encoder.hpp"
#include "../../decoder/CROM_decoder.hpp"
using namespace std;

int main() {
    double R = 0.01;
    int xdim = TEST_BLOCKLENGTH;
    double x_save[TEST_BLOCKLENGTH];
    bool verbose = false;

    double xhat[TEST_BLOCKLENGTH];

    CROM_encoder enc (xdim, R, verbose);
    string fname = "../input/x_input.txt";
    
    // read x
    enc.read_x(fname);
    // store original x
    enc.copy_x(x_save);

    cout << "Running CROM : Encoding" << endl;
    // encoding
    enc.run();

    int L = enc.get_L();
    int *m_array_copy = new int[L];
    enc.copy_m_array(m_array_copy);

    cout << "Running CROM : Decoding" << endl;
    // decoding
    int read_line_idx;
    for (read_line_idx=0; read_line_idx<TEST_BLOCKLENGTH; read_line_idx++) {
        xhat[read_line_idx] = 0;
    }
    CROM_decoder(xhat, xdim, L, m_array_copy, verbose);

    cout << "Comparing" << endl;
    double l2norm = 0;
    double diff;
    for (read_line_idx=0; read_line_idx<TEST_BLOCKLENGTH; read_line_idx++) {
        diff = x_save[read_line_idx] - xhat[read_line_idx];
        l2norm += (diff*diff);
    }
    double n = static_cast<double> (xdim);
    l2norm /= n;
    cout << "l2-norm at the end = " << l2norm << endl;
    delete[] m_array_copy;
    return 0;
}

