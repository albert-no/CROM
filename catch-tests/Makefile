all:
	g++ -std=c++11 test_crom_encoder.cpp test_fftw.cpp test_crom_big.cpp\
            test_matrix_multiplication.cpp test_crom_full.cpp\
            ../decoder/crom_decoder.cpp ../encoder/crom_encoder.cpp\
            ../utils/matrix_multiplication.cpp ../utils/crom_util.cpp\
            -L/usr/local/lib -I/usr/local/include -lfftw3 -o test_crom.o -lm
