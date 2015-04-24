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
  
  HEAD should be now safe enough to be used

  b) Download the armadillo library from http://arma.sourceforge.net/download.html 
  
  c) Modify config.mk to indicate the correct location of this armadillo library, if needed.
  
  d) Set the environment variable LD_LIBRARY_PATH to point to the armadillo library, id needed.

All codes are provided with some basic help, just use the 
--help option:

usage: ./readrootfile [options] rootfile 

 -h, --help               : display this help and exit
 -b, --bank-stubs         : extract bankstubs
 -n, --bank-stubs-new     : extract bankstubs new version
 -e, --not-use-charge     : if new version is true by default read also charge
 -l, --l1tk-stubs         : extract l1tkstubs
 -m, --max-tracks[=value] : max value of tracks to be extracted

split r-phi and r-z plane:

usage: ./generatepca_split [options] coordinatesfile 

 -h, --help                      : display this help and exit
 -v, --version                   : print version and exit
 -j, --jump-tracks               : generate the constants using only even tracks
 -p, --dump-allcoords            : dump all stub coordinates to a file
 -z, --rz-plane                  : use rz plane view
 -r, --rphi-plane                : use r-phi plane view
 -e, --not-use-charge            : do not read charge from coordinatesfile, by default 
                                   we will  use it if rphi-plane has been selected
 -g, --charge-sign=[+/-]         : use only + particle or - paricle 
 -t, --eta-range="etamin,etamax" : specify the eta range to use 
 -x, --exclude-s-module          : exclude S-module (last three layer) so 6 coordinates inseatd of 12 
 -d, --use-d0                    : use also d0 param in r-phi plane 

usage: ./fitpca_split [options] coordinatesfile 

 -h, --help                      : display this help and exit
 -V, --verbose                   : verbose option on
 -v, --version                   : print version and exit
 -c, --cmtx=[fillename]          : CMTX filename [default is c.[rz/rphi].bin]
 -q, --qvct=[fillename]          : QVCT filename [default is q.[rz/rphi].bin]
 -j, --jump-tracks               : perform the fittin only for odd tracks
 -z, --rz-plane                  : use rz plane view
 -r, --rphi-plane                : use r-phi plane view
 -e, --not-use-charge            : do not read charge from coordinatesfile, by default 
                                   and use it if rphi-plane has been selected
 -g, --charge-sign=[+/-]         : use only + particle or - paricle 
 -t, --eta-range="etamin;etamax" : specify the eta range to use 
 -x, --exclude-s-module          : exclude S-module (last three layer) so 6 coordinates inseatd of 12 
 -d, --use-d0                    : use also d0 param in r-phi plane 

old approach to use specific sub-tower and / or sub-sub-tower 

usage: ./bin/generatepca [options] coordinatesfile 

 -h, --help                 : display this help and exit
 -V, --verbose              : verbose option on
 -v, --version              : print version and exit
 -f, --fast                 : do not perfomr pca only diag matrix
 -i, --seg-id               : use SegId instead of coordinates
 -s, --bigger-sub-tower     : use values of the bigger sub-tower
                              (connot be used with bigger-subladder)
 -l, --bigger-sub-sub-tower : use values of the bigger sub-sub-tower 
                              (connot be used with bigger-subsector)
 -a, --all-sub-tower        : generate the constants for all subsectors
 -r, --all-sub-sub-tower    : generate the constants for all subladders
 -j, --jump-tracks          : generate the constants using only even tracks
 -p, --dump-allcoords       : dump all stub coordinates to a file

usage: ./bin/fitpca [options] coordinatesfile 

 -h, --help               : display this help and exit
 -V, --verbose            : verbose option on
 -v, --version            : print version and exit
 -c, --cmtx=[fillename]   : CMTX filename [default is c.bin]
 -q, --qvct=[fillename]   : QVCT filename [default is q.bin]
 -i, --seg-id             : use SegId instead of coordinates
 -s, --subsector=[subsec] : by default use values of the bigger subsector
                            with this option you can speficy to perform 
                            prediction for subsector subsec 
 -l, --subladder=[subld]  : by default use values of the bigger subladder 
                            with this option you can speficy to perform 
                            prediction for subladder subld 
 -a, --all-subsectors     : perform the fitting for all subsectors
 -r, --all-subladders     : perform the fitting for all subladders
 -j, --jump-tracks        : perform the fittin only for odd tracks


Example: 

TODO


