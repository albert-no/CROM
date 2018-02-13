all:
	g++ -std=c++11 test_cromq.cpp\
	    ../../decoder/crom_decoder.cpp ../../encoder/crom_encoder.cpp\
	    ../../utils/matrix_multiplication.cpp ../../utils/crom_util.cpp\
	    ../cromq_encoder.cpp ../cromq_decoder.cpp ../cromq_util.cpp\
	    -L/usr/local/lib -I/usr/local/include -I/home/noa/mycode/CROM/cromq -lfftw3 -lm -o test_CROMq.o
