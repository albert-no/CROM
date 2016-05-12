#include <fstream>
#define TEST_BLOCKLENGTH 16
#define TEST_HALFBLOCKLENGTH 16

#include "../encoder/CROM_encoder.hpp"
using namespace std;

int main() {
    double x[TEST_BLOCKLENGTH];
    double thetas[TEST_HALFBLOCKLENGTH];
    double R = 1;
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
    x_infile.open("x_input.txt");
    
    int read_line_idx;
    for (read_line_idx=0; read_line_idx<TEST_BLOCKLENGTH; read_line_idx++) {
        x_infile >> x[read_line_idx];
    }
    x_infile.close();

    cout << "Running CROM" << endl;
    CROM_encoder(x, xdim, L, m_array, true);

    int m_iter_idx;
    for (m_iter_idx=0; m_iter_idx<L; m_iter_idx++) {
        cout << m_array[m_iter_idx] << endl;
    }
    free(m_array);
    return 0;
}

