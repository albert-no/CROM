all:
	g++ -std=c++11 test_CROM_encoder.cpp test_fftw.cpp\
            test_matrix_multiplication.cpp test_CROM_full.cpp\
            ../decoder/CROM_decoder.cpp ../encoder/CROM_encoder.cpp\
            ../utils/matrix_multiplication.cpp ../utils/CROM_util.cpp\
            -L/usr/local/lib -I/usr/local/include -lfftw3 -o test_CROM.o -lm
