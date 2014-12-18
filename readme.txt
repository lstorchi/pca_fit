README (Loriano Storchi) loriano@storchi.org
---------------------------------------------------------------------
- basic pca fitting methods. Here  you  find  the  rootfile  reader 
  as well as the generatepca constants code and the code performing 
  the fitting.
- IMPORTANTE NOTE: This is just a test code for internal use only, 
  nothing more than a basic rapid quick and dirty prototyping. In case 
  we will never need a proper tool we will need to think about and 
  totally restructuring the code in a proper OO fashion. 
- In case you want to test it use: rel_0_2_0
 

All codes are provided with some basic help, just use the 
--help option,

usage: ./readrootfile [options] rootfile 

 -h, --help               : display this help and exit
 -b, --bank-stubs         : extract bankstubs
 -l, --l1tk-stubs         : extract l1tkstubs
 -m, --max-tracks[=value] : max value of tracks to be extracted

usage: ./generatepca [options] coordinatesfile 

 -h, --help               : display this help and exit
 -V, --verbose            : verbose option on
 -v, --version            : print version and exit
 -f, --fast               : do not perfomr pca only diag matrix
 -s, --bigger-subsector   : use values of the bigger subsector
                            (connot be used with bigger-subladder)
 -l, --bigger-subladder   : use values of the bigger subladder 
                            (connot be used with bigger-subsector)
 -a, --all-subsectors     : generate the constants for all subsectors
 -r, --all-subladders     : generate the constants for all subladders

usage: ./fitpca [options] coordinatesfile 

 -h, --help               : display this help and exit
 -V, --verbose            : verbose option on
 -v, --version            : print version and exit
 -c, --cmtx=[fillename]   : CMTX filename [default is c.bin]
 -q, --qvct=[fillename]   : QVCT filename [default is q.bin]
 -s, --subsector=[subsec] : by default use values of the bigger subsector
                            with this option you can speficy to perform 
                            prediction for subsector subsec 
 -l, --subladder=[subld]  : by default use values of the bigger subladder 
                            with this option you can speficy to perform 
                            prediction for subladder subld 
 -a, --all-subsectors     : perform the fitting for all subsectors
 -r, --all-subladders     : perform the fitting for all subladders
