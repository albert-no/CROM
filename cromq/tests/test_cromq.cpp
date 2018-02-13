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
    double R_dec = 0.04;
    double rd_param = 1.4;

    int id = 0;
    int xdim = TEST_BLOCKLENGTH;
    int num_x = 36;
    bool verbose = false;

    string name = "cromq_test";
    string fname = "sample.qscore";

    srand(SRAND_SEED);

    std::clock_t cromq_time;

    cromq_time = std::clock();

    // CROMq_encoder enc(name, fname, id, num_x, xdim, rd_param, R_enc, verbose);
    // enc.run();

    // cromq_time = std::clock() - cromq_time;

    // std::cout << "time = " << ((double)cromq_time / (double)CLOCKS_PER_SEC) << std::endl;

    // CROMq_decoder dec(name, fname, id, num_x, xdim, rd_param, R_dec, true);
    // dec.run();

    std::string ifname = "sample.qscore";
    std::string ofname = "cromq_test_id_0_out.txt";

    double distortion;
    distortion = compute_distortion(ifname, ofname, num_x, xdim);
    std::cout << "distortion = " << distortion << std::endl;
    return 0;
}

