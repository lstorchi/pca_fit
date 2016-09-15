CFLAGS +=-I../../include/ -std=c++11

all: coeff_file_gen translate_file

coeff_file_gen:
	g++ $(CFLAGS) -o ../../bin/coeff_file_gen coeff_file_gen.cc utility.cpp

translate_file:
	g++ $(CFLAGS) -o ../../bin/translate_file translate_file.cpp utility.cpp

clean:
	rm -f ../../bin/coeff_file_gen ../../bin/translate_file
