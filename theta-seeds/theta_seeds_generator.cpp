/*
   theta_seeds_generator.cpp
   Generate theta seeds
   Picked SEED_MAX as a prime number for no reason
*/
#include <fstream>
#include <cstdlib>
#include <cmath>

#define BLOCKLENGTH 65536
#define SEED_MAX 32573

using namespace std;

int main() {
    long write_line_idx;
    double rand_int;

    ofstream theta_seeds_file;
    theta_seeds_file.open("theta_seeds.txt");
    for (write_line_idx=0; write_line_idx<BLOCKLENGTH; write_line_idx++) {
        rand_int = rand() % SEED_MAX;
        // rand_int = write_line_idx;
        theta_seeds_file << rand_int << endl;
    }
    theta_seeds_file.close();

    return 0;
}
