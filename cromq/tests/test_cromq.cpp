#include <ctime>

#define TEST_BLOCKLENGTH 65536
#define TEST_HALFBLOCKLENGTH 32767
#define SRAND_SEED 5298
#define EPSILON 1e-6

#include "catch.hpp"
#include "../../encoder/crom_encoder.hpp"
#include "../../decoder/crom_decoder.hpp"
#include "../../utils/crom_util.hpp"
#include "../cromq_encoder.hpp"
//#include "../cromq_decoder.hpp"

using namespace std;

TEST_CASE("CROMq", "[CROM_q]") {
    double R = 1.0;
    double rd_param = 1.4;

    int xdim = TEST_BLOCKLENGTH;
    int num_x = 36;

    string name = "cromq_test";
    string fname = "sample.qscore";

    CROMq_encoder enc(name, fname, num_x, xdim, rd_param, R);


    SECTION( "compare final l2 norm after big run" ) {
        enc.run();
    }
}

