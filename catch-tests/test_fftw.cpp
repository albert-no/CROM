#include <fftw3.h>
#include <iostream>
#include <vector>

#define TEST_BLOCKLENGTH 8
#define EPSILON 1e-6

#include "catch.hpp"
#include "../utils/CROM_util.hpp"
#include "../utils/FastDctLee.hpp"

using namespace std;

TEST_CASE("fftw small test", "[fftw]") {
    fftw_plan p;
    // double x[TEST_BLOCKLENGTH] = {
    std::vector<double> x = {
        1.1, 1.2, -1.0, 1.0,
        1.0, 1.6, 1.0, 1.0};
    std::vector<double> xout(TEST_BLOCKLENGTH);
    int i;
    int xdim = TEST_BLOCKLENGTH;
    int x_iter;

    SECTION( "test 1" ) {
        // check output after FFT
        p = fftw_plan_r2r_1d(xdim, x.data(), xout.data(), FFTW_REDFT10, FFTW_ESTIMATE);
        fftw_execute(p);
        //double xout_expected[TEST_BLOCKLENGTH] = {
        std::vector<double> xout_expected = {
            13.8, -2.3602203107,
            1.4093628901, 5.1883412518,
            1.8384776311, -1.29566974,
            -2.8798778176, -4.5068520128};
        for (x_iter=0; x_iter<xdim; x_iter++) {
            CHECK( xout[x_iter] == Approx(xout_expected[x_iter]).epsilon(EPSILON) );
        }

        // check after normalize
        normalize_then_copy_vector(x, xout, xdim);
        //double x_norm_expected[TEST_BLOCKLENGTH] = {
        std::vector<double> x_norm_expected = {
            2.4395183951, -0.5900550777,
            0.3523407225, 1.2970853129,
            0.4596194078, -0.323917435,
            -0.7199694544, -1.1267130032};
        for (x_iter=0; x_iter<xdim; x_iter++) {
            CHECK( x[x_iter] == Approx(x_norm_expected[x_iter]).epsilon(EPSILON) );
        }
        fftw_destroy_plan(p);

        // check after unnormalize (before IFFT)
        unnormalize_vector(x, xdim);
        // double x_unnorm_expected[TEST_BLOCKLENGTH] = {
        std::vector<double> x_unnorm_expected = {
            0.8625, -0.1475137694,
            0.0880851806, 0.3242713282,
            0.1149048519, -0.0809793588,
            -0.1799923636, -0.2816782508};
        for (x_iter=0; x_iter<xdim; x_iter++) {
            CHECK( x[x_iter] == Approx(x_unnorm_expected[x_iter]).epsilon(EPSILON) );
        }

        // check after IFFT
        p = fftw_plan_r2r_1d(xdim, x.data(), xout.data(), FFTW_REDFT01, FFTW_ESTIMATE);
        fftw_execute(p);
        copy_vector(x, xout, xdim);
        double x_final_expected[TEST_BLOCKLENGTH] = {
        1.1, 1.2, -1.0, 1.0,
        1.0, 1.6, 1.0, 1.0};
        for (x_iter=0; x_iter<xdim; x_iter++) {
            CHECK( x[x_iter] == Approx(x_final_expected[x_iter]).epsilon(EPSILON) );
        }
    }
}
