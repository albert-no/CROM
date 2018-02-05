all:
	g++ -std=c++11 test_CROM_encoder.cpp\
            ../decoder/CROM_decoder.cpp ../encoder/CROM_encoder.cpp\
            ../utils/matrix_multiplication.cpp ../utils/CROM_util.cpp\
            ../utils/FastDctFft.cpp ../utils/FftRealPair.cpp\
            -L/usr/local/lib -I/usr/local/include -lfftw3 -o test_CROM.o -lm
