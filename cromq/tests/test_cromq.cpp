#include <ctime>

#define TEST_BLOCKLENGTH 65536
#define TEST_HALFBLOCKLENGTH 32768
#define SRAND_SEED 5298
#define EPSILON 1e-6

#include "catch.hpp"
#include "../../encoder/crom_encoder.hpp"
#include "../../decoder/crom_decoder.hpp"
#include "../../utils/crom_util.hpp"
#include "../cromq_encoder.hpp"
//#include "../cromq_decoder.hpp"

using namespace std;

int main() {
    double R = 1;
    double rd_param = 1.4;

    int xdim = TEST_BLOCKLENGTH;
    int num_x = 36;
    bool verbose = false;

    string name = "cromq_test";
    string fname = "sample.qscore";

    srand(SRAND_SEED);

    std::clock_t cromq_time;

    cromq_time = std::clock();

    CROMq_encoder enc(name, fname, num_x, xdim, rd_param, R, verbose);

    enc.run();

    cromq_time = std::clock() - cromq_time;
    std::cout << "time = " << cromq_time << std::endl;
    return 0;
}

