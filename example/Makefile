all:
	g++ -std=c++17 main.cpp\
	    ../decoder/crom_decoder.cpp ../encoder/crom_encoder.cpp\
	    ../utils/matrix_multiplication.cpp ../utils/crom_util.cpp\
	    ../cromq/cromq_encoder.cpp ../cromq/cromq_decoder.cpp ../cromq/cromq_util.cpp\
	    ../cromq/cromq_wrapper.cpp\
	    -L/usr/local/lib -I/usr/local/include -I ../cromq/ -lfftw3 -lm -o crom_example.o
