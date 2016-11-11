CFLAGS +=-I../../include/ -std=c++11

all: coeff_file_gen translate_file add_eta_pt_ranges merge_missed_bin

coeff_file_gen:
	g++ $(CFLAGS) -o ../../bin/coeff_file_gen coeff_file_gen.cc utility.cpp

merge_missed_bin:
	g++ $(CFLAGS) -o ../../bin/merge_missed_bin merge_missed_bin.cpp utility.cpp

translate_file:
	g++ $(CFLAGS) -o ../../bin/translate_file translate_file.cpp utility.cpp

add_eta_pt_ranges:
	g++ $(CFLAGS) -o ../../bin/add_eta_pt_ranges add_eta_pt_ranges.cpp utility.cpp

clean:
	rm -f ../../bin/coeff_file_gen ../../bin/translate_file ../../bin/add_eta_pt_ranges \
	  ../../bin/merge_missed_bin
