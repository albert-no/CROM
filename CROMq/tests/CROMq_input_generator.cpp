// CROMq_input_generator.cpp
#include <cmath>
#include <cstdlib>
#include <fstream>

#define row_dim 655360
#define col_dim 102

using namespace std;
int main() {
    int row_iter, col_iter, col_idx;
    double uni_rand1, uni_rand2;
    double norm_rand1, norm_rand2;
    double std[col_dim];

    uni_rand1 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    std[col_dim-1] = uni_rand1;

    for (col_iter=col_dim-2; col_iter>=0; col_iter--) {
        uni_rand1 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        std[col_iter] = std[col_iter+1] + uni_rand1;
    }
    ofstream stdfile;
    stdfile.open("CROMq_test_stds.txt");
    for (col_iter=0; col_iter<col_dim; col_iter++) {
        stdfile << std[col_iter] << endl;
    }
    stdfile.close();

    ofstream ofile;
    ofile.open("CROMq_test_input.txt");

    for (row_iter=0; row_iter<row_dim; row_iter++) {
        for (col_iter=0; col_iter<col_dim/2; col_iter++) {
            uni_rand1 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            uni_rand2 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            col_idx = col_iter * 2;
            norm_rand1 = std[col_idx] * sqrt(-2*log(uni_rand1)) * cos(2*M_PI*uni_rand2);
            col_idx += 1;
            norm_rand2 = std[col_idx] * sqrt(-2*log(uni_rand1)) * sin(2*M_PI*uni_rand2);
            ofile << norm_rand1 << " " << norm_rand2 << " ";
        }
        ofile << endl;
    }
    ofile.close();
    return 0;
}
