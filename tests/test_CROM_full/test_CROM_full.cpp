#include <fstream>
#define TEST_BLOCKLENGTH 65536
#define TEST_HALFBLOCKLENGTH 32768

#include "../../encoder/CROM_encoder.hpp"
#include "../../decoder/CROM_decoder.hpp"
using namespace std;

int main() {
    double x[TEST_BLOCKLENGTH];
    double x_save[TEST_BLOCKLENGTH];
    double xhat[TEST_BLOCKLENGTH];
    double thetas[TEST_HALFBLOCKLENGTH];
    double R = 2;
    double doubleL;
    int xdim = TEST_BLOCKLENGTH;
    double n = static_cast<double> (xdim);
    int L;

    doubleL = n*R / log(n);
    L = static_cast<int> (doubleL);
    cout << "creating m array" << endl;
    cout << "L = " << L << endl;
    int *m_array = (int*) malloc(sizeof(int)*(L+1));
    
    cout << "reading x input" << endl;
    ifstream x_infile;
    x_infile.open("../input/x_input.txt");
    
    // read x and initialize xhat
    int read_line_idx;
    for (read_line_idx=0; read_line_idx<TEST_BLOCKLENGTH; read_line_idx++) {
        x_infile >> x[read_line_idx];
        x_save[read_line_idx] = x[read_line_idx];
        xhat[read_line_idx] = 0;
    }

    cout << "Running CROM" << endl;
    bool verbose = false;
    cout << "Encoding" << endl;
    CROM_encoder(x, xdim, L, m_array, verbose);

    cout << "Decoding" << endl;
    CROM_decoder(xhat, xdim, L, m_array, verbose);

    cout << "Comparing" << endl;
    double l2norm = 0;
    double diff;
    for (read_line_idx=0; read_line_idx<TEST_BLOCKLENGTH; read_line_idx++) {
        diff = x_save[read_line_idx] - xhat[read_line_idx];
        l2norm += (diff*diff);
    }
    l2norm /= n;
    cout << "l2-norm at the end = " << l2norm << endl;
    free(m_array);
    return 0;
}

