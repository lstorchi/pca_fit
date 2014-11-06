README (Loriano Storchi) loriano@storchi.org
---------------------------------------------------------------------
- basic pca fitting methods. Here  you  find  the  rootfile  reader 
  as well as the generatepca constants code and the code performing 
  the fitting.
 

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


