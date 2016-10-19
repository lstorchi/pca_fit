#include "jobparams.h"

//tclap
#include <tclap/CmdLine.h>
using namespace TCLAP;

jobparams::jobparams(int argc, char** argv){
  


   try {
     // command line parser
     CmdLine cmd("Command option", ' ', "0.9");

     ValueArg<std::string> option("c","case","type of analysis (rates/sectors/rate_n_sec/sec_n_test/stub_eff/PR_eff)",
				false, "rates", "string");
     cmd.add(option);

     ValueArg<std::string> inputfile("i","input","path and name of the input file",
				false, "/scratch/viret/data.root", "string");
     cmd.add(inputfile);

     ValueArg<int> nsec("n","nsec","number of the sector to filter",
			false, 0, "int");
     cmd.add(nsec);

     ValueArg<int> lim("l","lim","the minimum number of layer/disks with stubs in the sector (never use less than 4 here)",
			false, 4, "int");
     cmd.add(lim);

     ValueArg<std::string> outfile("o","output","name of the output file",
				false, "output.root", "string");
     cmd.add(outfile);

     ValueArg<std::string> sector("s","sector","path and name of the sector file",
				  false, "/scratch/viret/data.csv", "string");
     cmd.add(sector);

     // parse
     cmd.parse(argc, argv);
     
     m_inputfile    = inputfile.getValue();
     m_outfile      = outfile.getValue();
     m_sector       = sector.getValue();
     m_opt          = option.getValue();
     m_nsec         = nsec.getValue();
     m_lim          = lim.getValue();
   }
   catch (ArgException &e){ // catch exception from parse
     std::cerr << "ERROR: " << e.error() << " for arg " << e.argId()  << std::endl;
     abort();
   }
}

