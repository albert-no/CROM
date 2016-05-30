#include <iostream>
#include <fftw3.h>
#include "../../utils/CROM_util.hpp"
#define TEST_BLOCKLENGTH 16
using namespace std;

int main() {
    fftw_plan p;
    double x[TEST_BLOCKLENGTH] = {1.1, 1.0, 1.0, 1.0,
                                  1.0, 1.0, 1.0, 1.0,
                                  1.0, 1.0, 1.0, 1.1,
                                  1.0, 1.0, 1.0, 1.0};
    double xout[TEST_BLOCKLENGTH];
    int i;
    int xdim = TEST_BLOCKLENGTH;
    print_vector(x, xdim);
    cout << "l2_norm of x = " << compute_l2(x, xdim) << endl;

    // test encoding (fft)
    p = fftw_plan_r2r_1d(xdim, x, xout, FFTW_REDFT10, FFTW_ESTIMATE);
    fftw_execute(p);
    cout << "FFT" << endl;
    print_vector(xout, xdim);
    // normalize
    normalize_then_copy_vector(x, xout, xdim);
    cout << "Normalization" << endl;
    print_vector(x, xdim);
    cout << "l2_norm of x = " << compute_l2(x, xdim) << endl;
    fftw_destroy_plan(p);

    // test decoding (ifft)
    cout << endl << "Unnormalization" << endl;
    // unnormalize
    unnormalize_vector(x, xdim);
    print_vector(x, xdim);

    p = fftw_plan_r2r_1d(xdim, x, xout, FFTW_REDFT01, FFTW_ESTIMATE);
    fftw_execute(p);
    cout << "IFFT" << endl;
    // normalize
    copy_vector(x, xout, xdim);
    print_vector(x, xdim);
    cout << "l2_norm of x = " << compute_l2(x, xdim) << endl;
    fftw_destroy_plan(p);
}
