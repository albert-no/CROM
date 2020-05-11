#include <ctime>
#include <vector>

#define TEST_BLOCKLENGTH 65536
#define TEST_HALFBLOCKLENGTH 32767
#define SRAND_SEED 5298
#define EPSILON 1e-6

#include "catch.hpp"
#include "../encoder/crom_encoder.hpp"
#include "../decoder/crom_decoder.hpp"
#include "../utils/crom_util.hpp"
using namespace std;

TEST_CASE("CROM_full R0.01", "[CROM_full]") {
    double R = 0.01;
    int xdim = TEST_BLOCKLENGTH;
    bool verbose = false;

    string name = "full_R0.01_test";

    std::vector<double> x(TEST_BLOCKLENGTH);
    std::vector<double> x_save(TEST_BLOCKLENGTH);
    std::vector<double> xhat(TEST_BLOCKLENGTH);

    srand(SRAND_SEED);

    CROM_encoder enc(name, xdim, R, verbose);
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
        CHECK( L == 40 );

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
        std::vector<double> l2_array_copy(L);
        enc.copy_l2_array(l2_array_copy);

        CHECK( l2dist == Approx(l2_array_copy[L-1]) );
        CHECK( l2dist == Approx(0.9868844746) );
    }
}

TEST_CASE("CROM_full R0.1", "[CROM_full]") {
    double R = 0.1;
    int xdim = TEST_BLOCKLENGTH;
    bool verbose = false;

    string name = "full_R0.1_test";

    std::vector<double> x(TEST_BLOCKLENGTH);
    std::vector<double> x_save(TEST_BLOCKLENGTH);
    std::vector<double> xhat(TEST_BLOCKLENGTH);

    srand(SRAND_SEED);

    CROM_encoder enc(name, xdim, R, verbose);
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
        CHECK( L == 409 );

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
        std::vector<double> l2_array_copy(L);
        enc.copy_l2_array(l2_array_copy);

        CHECK( l2dist == Approx(l2_array_copy[L-1]) );
        CHECK( l2dist == Approx(0.8963769396) );
    }
}

TEST_CASE("CROM_full R0.2", "[CROM_full]") {
    double R = 0.2;
    int xdim = TEST_BLOCKLENGTH;
    bool verbose = false;

    string name = "full_R0.2_test";

    std::vector<double> x(TEST_BLOCKLENGTH);
    std::vector<double> x_save(TEST_BLOCKLENGTH);
    std::vector<double> xhat(TEST_BLOCKLENGTH);

    srand(SRAND_SEED);

    CROM_encoder enc(name, xdim, R, verbose);
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
        CHECK( L == 819);

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
        std::vector<double> l2_array_copy(L);
        enc.copy_l2_array(l2_array_copy);

        CHECK( l2dist == Approx(l2_array_copy[L-1]) );
        CHECK( l2dist == Approx(0.8046017938) );
    }
}

TEST_CASE("CROM_full R0.5", "[CROM_full]") {
    double R = 0.5;
    int xdim = TEST_BLOCKLENGTH;
    bool verbose = false;

    string name = "full_R0.5_test";

    std::vector<double> x(TEST_BLOCKLENGTH);
    std::vector<double> x_save(TEST_BLOCKLENGTH);
    std::vector<double> xhat(TEST_BLOCKLENGTH);

    srand(SRAND_SEED);

    CROM_encoder enc(name, xdim, R, verbose);
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
        CHECK( L == 2048);

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
        std::vector<double> l2_array_copy(L);
        enc.copy_l2_array(l2_array_copy);

        CHECK( l2dist == Approx(l2_array_copy[L-1]) );
        CHECK( l2dist == Approx(0.5727847666) );
    }
}

TEST_CASE("CROM_full big", "[CROM_full]") {
    double R = 1.0;
    int xdim = TEST_BLOCKLENGTH;
    bool verbose = false;

    string name = "full_big_test";

    std::vector<double> x(TEST_BLOCKLENGTH);
    std::vector<double> x_save(TEST_BLOCKLENGTH);
    std::vector<double> xhat(TEST_BLOCKLENGTH);

    srand(SRAND_SEED);

    CROM_encoder enc(name, xdim, R, verbose);
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

    std::vector<double> l2_array_copy;

    double L01, L1, L2, L5;

    SECTION( "compare final l2 norm after big run" ) {
        // read and store original x
        enc.set_x(x);
        enc.copy_x(x_save);

        // check L
        int L = enc.get_L();
        l2_array_copy.resize(L);
        CHECK( L == 4096);

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
        enc.copy_l2_array(l2_array_copy);

        CHECK( l2dist == Approx(l2_array_copy[L-1]) );
        CHECK( l2dist == Approx(0.3240603953) );

        L01 = l2_array_copy[39];
        L1 = l2_array_copy[408];
        L2 = l2_array_copy[818];
        L5 = l2_array_copy[2047];
    }

    SECTION( "Partial Decoding R0.01" ) {
        // Rdec=0.01
        // decoding and set m_array
        int L = 40;
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

        printf("Checking l2dist of partial decoding R0.01\n");
        CHECK( l2dist == Approx(L01) );
        CHECK( l2dist == Approx(0.9868844746) );
    }

    SECTION( "Partial Decoding R0.1" ) {

        // Rdec=0.1
        // decoding and set m_array
        int L = 409;
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

        printf("Checking l2dist of partial decoding R0.1\n");
        CHECK( l2dist == Approx(L1) );
        CHECK( l2dist == Approx(0.8963769396) );
    }

    SECTION( "Partial Decoding R0.2" ) {

        // Rdec=0.2
        // decoding and set m_array
        int L = 819;
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

        printf("Checking l2dist of partial decoding R0.2\n");
        CHECK( l2dist == Approx(L2) );
        CHECK( l2dist == Approx(0.804601738) );
    }

    SECTION( "Partial Decoding R0.5" ) {

        // Rdec=0.5
        // decoding and set m_array
        int L = 2048;
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

        printf("Checking l2dist of partial decoding R0.5\n");
        CHECK( l2dist == Approx(L5) );
        CHECK( l2dist == Approx(0.5727847666) );
    }
}
