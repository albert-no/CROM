// test_cromq.test
#include <ctime>

#define TEST_BLOCKLENGTH 256
#define TEST_HALFBLOCKLENGTH 128
#define SRAND_SEED 5298
#define EPSILON 1e-5

#define CATCH_CONFIG_MAIN
#include "../../catch-tests/catch.hpp"
#include "../../encoder/crom_encoder.hpp"
#include "../../decoder/crom_decoder.hpp"
#include "../../utils/crom_util.hpp"
#include "../cromq_encoder.hpp"
#include "../cromq_decoder.hpp"
#include "../cromq_util.hpp"

using namespace std;

TEST_CASE("CROMq_full", "[CROMq_full]") {
    double R_enc = 5;
    double rd_param = 1.4;

    int num_dec_points = 5;
    int dec_idx;
    std::vector<double> R_dec = {1, 2, 3, 4, 5};

    int id = 0;
    int xdim = TEST_BLOCKLENGTH;
    int num_x = 36;
    bool verbose = false;

    string name = "cromq_test";
    string fname = "sample.qscore";
    string ofname;
    std::clock_t cromq_time;

    srand(SRAND_SEED);

    SECTION( "CROMq encoder test" ) {
        cromq_time = std::clock();

        CROMq_encoder enc(name, fname, id, num_x, xdim, rd_param, R_enc, verbose);
        enc.run();

        cromq_time = std::clock() - cromq_time;

        std::cout << "encoding time = " << ((double)cromq_time / (double)CLOCKS_PER_SEC) << std::endl;
    }

    SECTION( "CROMq decoder test") {
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
        }


        std::vector<double> D_array_expected = {4.65788, 2.68056, 2.19531, 2.04742, 2.00499};
        for (dec_idx=0; dec_idx<num_dec_points; dec_idx++) {
            CHECK( D_array[dec_idx] == Approx(D_array_expected[dec_idx]).epsilon(EPSILON) );
        }
    }
}

