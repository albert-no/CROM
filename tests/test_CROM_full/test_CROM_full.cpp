#define TEST_BLOCKLENGTH 65536
#define TEST_HALFBLOCKLENGTH 32768

#include "../../encoder/CROM_encoder.hpp"
#include "../../decoder/CROM_decoder.hpp"
#include "../../utils/CROM_util.hpp"
using namespace std;

int main() {
    double R = 0.1;
    int xdim = TEST_BLOCKLENGTH;
    double *x = new double[TEST_BLOCKLENGTH];
    double x_save[TEST_BLOCKLENGTH];
    bool verbose = false;

    double xhat[TEST_BLOCKLENGTH];
    string name = "full_test";

    cout << "Running CROM : Encoding" << endl;

    CROM_encoder enc (name, xdim, R, verbose);
    string string_fname = "../input/x_input.txt";

    int line_idx;
    ifstream x_infile;
    x_infile.open(string_fname.c_str());
    for (line_idx=0; line_idx<xdim; line_idx++) {
        x_infile >> x[line_idx];
    }
    x_infile.close();
    
    // read x
    enc.set_x(x);
    delete[] x;
    // store original x
    enc.copy_x(x_save);

    // encoding
    enc.run();

    int L = enc.get_L();
    int *m_array_copy = new int[L];
    // read m_array
    enc.copy_m_array(m_array_copy);

    cout << "Running CROM : Decoding" << endl;
    CROM_decoder dec (name, xdim, L, verbose);

    // set m_array
    dec.set_m_array(m_array_copy);

    // decoding
    dec.run();
    dec.copy_x_hat(xhat);

    // write x_hat into txt file
    dec.write_x_hat();

    cout << "Comparing" << endl;
    double l2dist;
    l2dist = compute_l2_dist(x_save, xhat, xdim);
    l2dist /= static_cast<double> (xdim);
    cout << "l2-dist at the end = " << l2dist << endl;
    delete[] m_array_copy;
    return 0;
}

