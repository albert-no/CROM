/*
   test_input_generator.cpp
   Generate input vector X and input thetas for test
*/
#include <fstream>
#include <cstdlib>
#include <cmath>

#define BLOCKLENGTH 65536
#define HALFBLOCKLENGTH 32768

using namespace std;

int main() {
    long write_line_idx;
    double uni_rand;

    ofstream x_outfile;
    x_outfile.open("x_input.txt");
    for (write_line_idx=0; write_line_idx<BLOCKLENGTH; write_line_idx++) {
        uni_rand = static_cast <float> (rand()) / static_cast <float> (RAND_MAX); 
        uni_rand -= 0.5;
        uni_rand *= sqrt(12);
        x_outfile << uni_rand << endl;
    }
    x_outfile.close();

    ofstream theta_outfile;
    theta_outfile.open("theta_input.txt");
    for (write_line_idx=0; write_line_idx<HALFBLOCKLENGTH; write_line_idx++) {
        uni_rand = static_cast <float> (rand()) / static_cast <float> (RAND_MAX); 
        theta_outfile << uni_rand*M_PI_2 << endl;
    }
    theta_outfile.close();
    return 0;
}
