/*
   test_matrix_multiplication.cpp
*/
#include "../utils/matrix_multiplication.hpp"
#include <fstream>

#define BLOCKLENGTH 65536
#define HALFBLOCKLENGTH 32768

using namespace std;

int main() {
    double x[BLOCKLENGTH];
    double thetas[HALFBLOCKLENGTH];

    ifstream x_infile;
    x_infile.open("x_input.txt");

    long read_line_idx;
    for (read_line_idx=0; read_line_idx<BLOCKLENGTH; read_line_idx++) {
        x_infile >> x[read_line_idx];
    }
    x_infile.close();

    ifstream theta_infile;
    theta_infile.open("theta_input.txt");
    for (read_line_idx=0; read_line_idx<HALFBLOCKLENGTH; read_line_idx++) {
        theta_infile >> thetas[read_line_idx];
    }
    theta_infile.close();
    
    butterfly_matrix_multiplication(x,
                                    thetas,
                                    HALFBLOCKLENGTH,
                                    0,
                                    BLOCKLENGTH,
                                    0,
                                    HALFBLOCKLENGTH,
                                    0);
    
    int iter_idx;
    ofstream x_outfile;
    x_outfile.open("x_output.txt");
    for (iter_idx=0; iter_idx<BLOCKLENGTH; iter_idx++) {
        x_outfile << x[iter_idx] << endl;
    }
    x_outfile.close();
    return 0;
}
