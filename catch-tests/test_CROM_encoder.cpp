// test_CROM_encoder.cpp
#include <limits>
#include <vector>
#define EPSILON 1e-6

#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../encoder/CROM_encoder.hpp"

TEST_CASE("small test", "[CROM_encoder]") {
    int x_dim = 8;
    std::string name = "basic_test";
    double R = 1;
    bool verbose = true;

    int x_iter; 
    CROM_encoder enc (name, x_dim, R, verbose);

    SECTION( "test 1" ) {
        // test set_x and copy_x
        double x[8] = {2.0, 1.2, 3.2, 1.1, 0.4, 0.7, 5.0, 10.2};
        enc.set_x(x);
        
        double x_copy[8];
        enc.copy_x(x_copy);

        for (x_iter=0; x_iter<x_dim; x_iter++) {
            CHECK( x[x_iter] == Approx(x_copy[x_iter]).epsilon(EPSILON) );
        }

        // check get_L
        int L = enc.get_L();
        CHECK( L == 3 );

        // run encoder
        enc.run();

        // check message array
        int *m_arr = new int[L];
        enc.copy_m_array(m_arr);
        std::vector<int> v_m_arr(m_arr, m_arr + L);
        std::vector<int> v_m_arr_expected {3, 3, 5};
        CHECK( v_m_arr == v_m_arr_expected );
        delete[] m_arr;

        // check x_remainder
        double x_rem[8] = {
            -2.901471318,
            -5.5536843625,
            -1.8770677173,
            -5.429310964,
            2.1008209927,
            2.5955557863,
            -1.9996325797,
            -2.6661145525};

        enc.copy_x(x_copy);
        for (x_iter=0; x_iter<x_dim; x_iter++) {
            CHECK( x_rem[x_iter] == Approx(x_copy[x_iter]).epsilon(EPSILON) );
        }
    }
}
