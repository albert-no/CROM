#define TEST_BLOCKLENGTH 65536
#define TEST_HALFBLOCKLENGTH 32768

#include "../../encoder/CROM_encoder.hpp"
using namespace std;

int main() {
    double R = 0.1;
    int xdim = TEST_BLOCKLENGTH;
    bool verbose = false;
    
    CROM_encoder enc (xdim, R, verbose);
    cout << "reading x input" << endl;
    string fname = "../input/x_input.txt";
    enc.read_x(fname);

    cout << "Running CROM" << endl;
    enc.run();

    enc.print_m_array();
    return 0;
}
