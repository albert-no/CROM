#include <iostream>
#include <fftw3.h>
#include "../../utils/CROM_util.hpp"
using namespace std;

int main() {
    fftw_plan p;
    double x[16] = {1.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    int i;
    for (i=0;i<16;i++) {
        cout << x[i] << " ";
    }
    cout << endl;
    cout << compute_l2(x, 16) << endl;
    p = fftw_plan_r2r_1d(16, x, x, FFTW_REDFT10, FFTW_ESTIMATE);
    fftw_execute(p);
    normalize_vector(x, 16);
    for (i=0;i<16;i++) {
        cout << x[i] << " ";
    }
    cout << endl;
    cout << compute_l2(x, 16) << endl;
    fftw_destroy_plan(p);
}

