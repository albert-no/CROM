#include <ctime>

#define TEST_BLOCKLENGTH 65536
#define TEST_HALFBLOCKLENGTH 32768
#define SRAND_SEED 5298
#define EPSILON 1e-6

#include "../../encoder/crom_encoder.hpp"
#include "../../decoder/crom_decoder.hpp"
#include "../../utils/crom_util.hpp"
#include "../cromq_encoder.hpp"
#include "../cromq_decoder.hpp"
#include "../cromq_util.hpp"

using namespace std;

int main() {
    double R_enc = 0.05;
    double rd_param = 1.4;

    int num_dec_points = 5;
    int dec_idx;
    std::vector<double> R_dec = {0.01, 0.02, 0.03, 0.04, 0.05};

    int id = 0;
    int xdim = TEST_BLOCKLENGTH;
    int num_x = 36;
    bool verbose = false;

    string name = "cromq_test";
    string fname = "sample.qscore";
    string ofname;

    srand(SRAND_SEED);

    std::clock_t cromq_time;

    cromq_time = std::clock();

    CROMq_encoder enc(name, fname, id, num_x, xdim, rd_param, R_enc, verbose);
    enc.run();

    cromq_time = std::clock() - cromq_time;

    std::cout << "encoding time = " << ((double)cromq_time / (double)CLOCKS_PER_SEC) << std::endl;

    std::vector<double> D_array(num_dec_points);
    for (dec_idx=0; dec_idx<num_dec_points; dec_idx++) {
        cromq_time = std::clock();
        CROMq_decoder dec(name, fname, id, num_x, xdim, rd_param, R_dec[dec_idx], false);
        dec.run();

        cromq_time = std::clock() - cromq_time;
        std::cout << "decoding time = " << ((double)cromq_time / (double)CLOCKS_PER_SEC) << std::endl;

        ofname = get_ofname(name, id, R_dec[dec_idx]);
        D_array[dec_idx] = compute_distortion(fname, ofname, num_x, xdim);
        std::cout << "R_dec = " << R_dec[dec_idx] << ", distortion = " << D_array[dec_idx] << std::endl;
        // 11.9941
        // 11.9209
        // 11.8412
        // 11.7614
        // 11.6736
    }

    return 0;
}

