# include <vector>

#define BLOCKLENGTH 16
#define HALFBLOCKLENGTH 8
#define EPSILON 1e-6

#include "catch.hpp"
#include "../utils/matrix_multiplication.hpp"

using namespace std;

TEST_CASE("matrix multiplication small test", "[matrix_multiplication]") {
    double x[BLOCKLENGTH] = {
        1.1, 1.0, 1.0, -2.1,
        0.1, 4.0, -1.0, 4.5,
        1.0, 2.0, 1.3, 0.1,
        -0.2, -0.3, 0.4, 2.2};
    double x_final_expected[BLOCKLENGTH] = {
        1.1, 1.0, 1.0, -2.1,
        0.1, 4.0, -1.0, 4.5,
        1.0, 2.0, 1.3, 0.1,
        -0.2, -0.3, 0.4, 2.2};

    long x_iter;

    std::vector<double> thetas = {
        0.1, 0.2, 0.3, 0.4,
        0.5, -0.2, -0.7, 2.0};

    std::vector<double> thetas_inv = {
        -0.1, -0.2, -0.3, -0.4,
        -0.5, 0.2, 0.7, -2.0};

    SECTION( "mat_idx = 0" ) {
        int mat_idx = 0;

        // matrix multiplication
        double x_out_expected[BLOCKLENGTH] = {
            0.9946711652, 0.5827279163,
            0.5711602205, -1.9731699216,
            0.1836433639, 3.8606655121,
            -0.5071551124, -3.8731151035,
            1.1048209236, 2.1588024865,
            1.5374576425, -0.7256724194,
            -0.1275739585, -1.0886972965,
            0.9501545622, 3.1763153803};

        butterfly_matrix_multiplication(x,
                                        thetas,
                                        HALFBLOCKLENGTH,
                                        0,
                                        0,
                                        mat_idx);
        for (x_iter=0; x_iter<BLOCKLENGTH; x_iter++) {
            CHECK( x[x_iter] == Approx(x_out_expected[x_iter]).epsilon(EPSILON) );
        }
        
        // inverse matrix multiplication
        butterfly_matrix_multiplication(x,
                                        thetas_inv,
                                        HALFBLOCKLENGTH,
                                        0,
                                        0,
                                        mat_idx);
        for (x_iter=0; x_iter<BLOCKLENGTH; x_iter++) {
            CHECK( x[x_iter] == Approx(x_final_expected[x_iter]).epsilon(EPSILON) );
        }
    }

    SECTION( "mat_idx = 2" ) {
        int mat_idx = 2;

        // matrix multiplication
        double x_out_expected[BLOCKLENGTH] = {
            0.9946711652, 1.3972721725,
            1.1048209236, -1.8594704827,
            0.3910538556, 1.9318614356,
            -0.9257844685, 5.7024478422,
            0.2543293617, 1.9800000888,
            1.6202828691, -0.2993320038,
            0.1047186374, -1.8756102881,
            0.4347804124, -1.1883122685};
        butterfly_matrix_multiplication(x,
                                        thetas,
                                        HALFBLOCKLENGTH,
                                        0,
                                        0,
                                        mat_idx);
        for (x_iter=0; x_iter<BLOCKLENGTH; x_iter++) {
            CHECK( x[x_iter] == Approx(x_out_expected[x_iter]).epsilon(EPSILON) );
        }
        
        // inverse matrix multiplication
        butterfly_matrix_multiplication(x,
                                        thetas_inv,
                                        HALFBLOCKLENGTH,
                                        0,
                                        0,
                                        mat_idx);
        for (x_iter=0; x_iter<BLOCKLENGTH; x_iter++) {
            CHECK( x[x_iter] == Approx(x_final_expected[x_iter]).epsilon(EPSILON) );
        }
    }
}
