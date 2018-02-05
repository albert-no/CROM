// test_CROM_encoder.cpp
#include <limits>
#include <vector>
#define EPSILON 1e-6

#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../encoder/CROM_encoder.hpp"

TEST_CASE("CROM_encoder small test", "[CROM_encoder]") {
    int x_dim = 8;
    std::string name = "basic_test";
    double R = 1;
    bool verbose = true;

    int x_iter; 
    CROM_encoder enc (name, x_dim, R, verbose);

    SECTION( "test 1" ) {
        // test set_x and copy_x
        std::vector<double> x = {
            1.17844, -0.365868,
            0.980685, 1.03383,
            1.42599, -1.04771,
            -0.570805, 0.929175};
        enc.set_x(x);
        
        std::vector<double> x_copy(x_dim);
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
        std::vector<int> v_m_arr(L);
        enc.copy_m_array(v_m_arr);
        std::vector<int> v_m_arr_expected {2, 7, 0};
        CHECK( v_m_arr == v_m_arr_expected );

        // check x_remainder
        std::vector<double> x_rem = {
            -0.4673188324,
            -1.3669546712,
            -0.9557097795,
            -0.1965084036,
            0.5728391238,
            0.2799289636,
            0.3629543445,
            -0.6427417279};
        enc.copy_x(x_copy);
        for (x_iter=0; x_iter<x_dim; x_iter++) {
            CHECK( x_rem[x_iter] == Approx(x_copy[x_iter]).epsilon(EPSILON) );
        }
    }
}
