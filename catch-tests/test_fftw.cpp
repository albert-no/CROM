#include <fftw3.h>
#include <iostream>
#include <vector>

#define TEST_BLOCKLENGTH 8
#define EPSILON 1e-6

#include "catch.hpp"
#include "../utils/CROM_util.hpp"
#include "../utils/FastDctFft.hpp"
using namespace std;

TEST_CASE("fftw small test", "[fftw]") {
    std::vector<double> x = {1.1, 1.2, -1.0, 1.0,
                             1.0, 1.6, 1.0, 1.0};
    int i;
    int xdim = TEST_BLOCKLENGTH;
    int x_iter;

    SECTION( "test 1" ) {
        // check output after FFT
        FastDctFft::transform(x);
        std::vector<double> xout_expected = {
            6.9, -1.1801101554,
            0.704681445, 2.5941706259,
            0.9192388155, -0.64783487,
            -1.4399389088, -2.2534260064};
        for (x_iter=0; x_iter<xdim; x_iter++) {
            CHECK( x[x_iter] == Approx(xout_expected[x_iter]).epsilon(EPSILON) );
        }

        // check after normalize
        normalize_vector(x, xdim);
        std::vector<double> x_norm_expected = {
            2.4395183951, -0.5900550777,
            0.3523407225, 1.2970853129,
            0.4596194078, -0.323917435,
            -0.7199694544, -1.1267130032};
        for (x_iter=0; x_iter<xdim; x_iter++) {
            CHECK( x[x_iter] == Approx(x_norm_expected[x_iter]).epsilon(EPSILON) );
        }

        // check after unnormalize (before IFFT)
        unnormalize_vector(x, xdim);
        std::vector<double> x_unnorm_expected = {
            1.725, -0.2950275388,
            0.1761703613, 0.6485426565,
            0.2298097039, -0.1619587175,
            -0.3599847272, -0.5633565016};
        for (x_iter=0; x_iter<xdim; x_iter++) {
            CHECK( x[x_iter] == Approx(x_unnorm_expected[x_iter]).epsilon(EPSILON) );
        }

        // check after IFFT
        FastDctFft::inverseTransform(x);
        double x_final_expected[TEST_BLOCKLENGTH] = {
            1.1, 1.2, -1.0, 1.0,
            1.0, 1.6, 1.0, 1.0};
        for (x_iter=0; x_iter<xdim; x_iter++) {
            CHECK( x[x_iter] == Approx(x_final_expected[x_iter]).epsilon(EPSILON) );
        }
    }
}
