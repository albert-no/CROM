#include <ctime>

#define TEST_BLOCKLENGTH 65536
#define TEST_HALFBLOCKLENGTH 32767
#define SRAND_SEED 5298
#define EPSILON 1e-6

#include "catch.hpp"
#include "../encoder/CROM_encoder.hpp"
#include "../decoder/CROM_decoder.hpp"
#include "../utils/CROM_util.hpp"
using namespace std;

TEST_CASE("CROM_full small big", "[CROM_full]") {
    double R = 1.0;
    int xdim = TEST_BLOCKLENGTH;
    srand(SRAND_SEED);
    double x[TEST_BLOCKLENGTH];
    double x_save[TEST_BLOCKLENGTH];
    bool verbose = false;
    double xhat[TEST_BLOCKLENGTH];
    string name = "full_big_test";

    CROM_encoder enc (name, xdim, R, verbose);
    int x_iter;
    double uniform1, uniform2;
    double normal1, normal2;
    for (x_iter=0; x_iter<xdim/2; x_iter++) {
        uniform1 = static_cast <double> (rand()) /static_cast <double> (RAND_MAX);
        uniform2 = static_cast <double> (rand()) /static_cast <double> (RAND_MAX);
        normal1 = sqrt(-2*log(uniform1)) * cos(2*M_PI*uniform2);
        normal2 = sqrt(-2*log(uniform1)) * sin(2*M_PI*uniform2);
        x[2*x_iter] = x_save[2*x_iter] = normal1;
        x[2*x_iter+1] = x_save[2*x_iter+1] = normal2;
    }


    SECTION( "compare final l2 norm after big run" ) {
        // read and store original x
        enc.set_x(x);
        enc.copy_x(x_save);

        // check L
        int L = enc.get_L();
        CHECK( L == 5909 );

        std::clock_t enc_time;
        printf("BIG test: xdim = %d, and R = %f\n", xdim, R);
        // encoding
        enc_time = std::clock();
        enc.run();
        enc_time = std::clock() - enc_time;
        printf("Encoding took %f seconds \n", ((double)enc_time)/(double)CLOCKS_PER_SEC);

        // write_m_array in binary file
        enc.write_m_array(true);

        // decoding and set m_array
        CROM_decoder dec (name, xdim, L, verbose);
        dec.read_m_array(true);

        std::clock_t dec_time;
        // decoding
        dec_time = clock();
        dec.run();
        dec_time = clock() - dec_time;
        printf("Decoding took %f seconds \n", ((double)dec_time)/(double)CLOCKS_PER_SEC);
        dec.copy_x_hat(xhat);

        // check write_x_hat
        dec.write_x_hat();

        // check final l2 distance between xhat and x
        double l2dist;
        l2dist = compute_l2_dist(x_save, xhat, xdim);
        l2dist /= static_cast<double> (xdim);
        double *l2_array_copy = new double[L];
        enc.copy_l2_array(l2_array_copy);

        CHECK( l2dist == Approx(l2_array_copy[L-1]) );
        delete[] l2_array_copy;
    }
}

