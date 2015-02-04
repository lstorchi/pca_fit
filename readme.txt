README (Loriano Storchi) loriano@storchi.org
---------------------------------------------------------------------
- basic pca fitting methods. Here  you  find  the  rootfile  reader 
  as well as the generatepca constants code and the code performing 
  the fitting.
- IMPORTANTE NOTE: This is just a test code for internal use only, 
  nothing more than a basic rapid quick and dirty prototyping (i.e.
  the code should be totally restructured in a proper OO fashion)
- also HEAD rev should be safe enough to be used. Step by step:

  a) dowload the code  

  $ git clone https://bitbucket.org/lstorchi/gf_fit PCA
  $ cd PCA 
  
  to use rel_0_2_0 

  $ git checkout rel_0_2_0       

  b) Download the armadillo library from http://arma.sourceforge.net/download.html 
  
  c) Modify config.mk to indicate the correct location of this armadillo library, if needed.
  
  d) Set the environment variable LD_LIBRARY_PATH to point to the armadillo library, id needed.

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
