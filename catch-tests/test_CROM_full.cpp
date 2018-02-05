#define TEST_BLOCKLENGTH 8
#define TEST_HALFBLOCKLENGTH 4
#define EPSILON 1e-6

#include "catch.hpp"
#include "../encoder/CROM_encoder.hpp"
#include "../decoder/CROM_decoder.hpp"
#include "../utils/CROM_util.hpp"
using namespace std;

TEST_CASE("CROM_full small test", "[CROM_full]") {
    double R = 1;
    int xdim = TEST_BLOCKLENGTH;
    bool verbose = false;
    string name = "full_test";
    std::vector<double> x = {
        -0.769811, 0.186957,
        -0.0782989, 0.446422,
        -0.4684, 0.0454221,
        1.56657, 1.44174};
    std::vector<double> x_save = {
        -0.769811, 0.186957,
        -0.0782989, 0.446422,
        -0.4684, 0.0454221,
        1.56657, 1.44174};
    std::vector<double> xhat_expected = {
        -0.4883347668, 0.3162271728,
        -0.0086511553, 0.375842363,
        -0.5104244858, 0.1219668209,
        0.23900423, 0.825249401};

    std::vector<int> v_m_arr_expected {3, 1, 1};
    std::vector<double> xhat(TEST_BLOCKLENGTH);

    CROM_encoder enc (name, xdim, R, verbose);
    int x_iter;
    
    SECTION( "basic test with copying m_array" ) {
        // read and store original x
        enc.set_x(x);
        enc.copy_x(x_save);

        // check L
        int L = enc.get_L();
        CHECK( L == 3 );

        // encoding and check m_array
        enc.run();
        std::vector<int> v_m_arr(L);
        enc.copy_m_array(v_m_arr);
        CHECK( v_m_arr == v_m_arr_expected );

        // check print_m_array
        printf("Test print_m_array\n");
        enc.print_m_array();

        // decoding and set m_array
        CROM_decoder dec (name, xdim, L, verbose);
        dec.set_m_array(v_m_arr);

        // decoding
        dec.run();
        dec.copy_x_hat(xhat);
        for (x_iter=0; x_iter< xdim; x_iter++) {
            CHECK( xhat[x_iter] == Approx(xhat_expected[x_iter]).epsilon(EPSILON) );
        }

        // check write_x_hat
        dec.write_x_hat();

        // check final l2 distance between xhat and x
        double l2dist;
        l2dist = compute_l2_dist(x_save, xhat, xdim);
        l2dist /= static_cast<double> (xdim);
        CHECK( l2dist == Approx(0.2819860781) );
    }

    SECTION( "write and read m_array using txt file" ) {
        // read and store original x
        enc.set_x(x);
        enc.copy_x(x_save);
        int L = enc.get_L();

        // encoding
        enc.run();

        // check write_m_array
        enc.write_m_array(false);  // txt file

        // decoding and set m_array
        CROM_decoder dec (name, xdim, L, verbose);
        dec.read_m_array(false);  // txt file

        // decoding
        dec.run();
        dec.copy_x_hat(xhat);
        for (x_iter=0; x_iter< xdim; x_iter++) {
            CHECK( xhat[x_iter] == Approx(xhat_expected[x_iter]).epsilon(EPSILON) );
        }

        // check write_x_hat
        dec.write_x_hat();

        // check final l2 distance between xhat and x
        double l2dist;
        l2dist = compute_l2_dist(x_save, xhat, xdim);
        l2dist /= static_cast<double> (xdim);
        CHECK( l2dist == Approx(0.2819860781) );
    }

    SECTION( "write and read m_array using binary file" ) {
        // read and store original x
        enc.set_x(x);
        enc.copy_x(x_save);
        int L = enc.get_L();

        // encoding and check m_array
        enc.run();

        // write_m_array in binary file
        enc.write_m_array(true);

        // decoding and set m_array
        CROM_decoder dec (name, xdim, L, verbose);
        dec.read_m_array(true);

        // decoding
        dec.run();
        dec.copy_x_hat(xhat);
        for (x_iter=0; x_iter< xdim; x_iter++) {
            CHECK( xhat[x_iter] == Approx(xhat_expected[x_iter]).epsilon(EPSILON) );
        }

        // check write_x_hat
        dec.write_x_hat();

        // check final l2 distance between xhat and x
        double l2dist;
        l2dist = compute_l2_dist(x_save, xhat, xdim);
        l2dist /= static_cast<double> (xdim);
        CHECK( l2dist == Approx(0.2819860781) );

        // compare final l2 distance and encoder's l2_array
        std::vector<double> l2_array_copy(L);
        enc.copy_l2_array(l2_array_copy);
        CHECK( l2dist == Approx(l2_array_copy[L-1]) );
    }
}

