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
  
  to use rel_0_2_0 (HEAD should be now safe enough to be used)  

  $ git checkout rel_0_2_0       

  b) Download the armadillo library from http://arma.sourceforge.net/download.html 
  
  c) Modify config.mk to indicate the correct location of this armadillo library, if needed.
  
  d) Set the environment variable LD_LIBRARY_PATH to point to the armadillo library, id needed.

All codes are provided with some basic help, just use the 
--help option:

usage: ./bin/readrootfile [options] rootfile 

 -h, --help               : display this help and exit
 -b, --bank-stubs         : extract bankstubs
 -n, --bank-stubs-new     : extract bankstubs new version
 -l, --l1tk-stubs         : extract l1tkstubs
 -m, --max-tracks[=value] : max value of tracks to be extracted

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

$ ./readrootfile -n DistrD0.root

$ ./generatepca -l ./bankstub.txt 
file has 57513 line 
file has 7189 entries 
Reading data from ./bankstub.txt file 
Writing parameters to files
Looking for bigger subladder
We  found 323 subladders 
Selected subladder 050531-060727-071027-081412-091812-102212 numevt: 566
Printout selected coordinates 
Writing extracted parameters to files
Perform PCA 
566 18
Writing Scoreplot
Eigenvalues: 
1 ==> 42.6131% value: 9.25344
2 ==> 35.3698% value: 7.68055
3 ==> 11.0057% value: 2.3899
4 ==> 8.00011% value: 1.73722
5 ==> 2.01361% value: 0.437254
6 ==> 0.611055% value: 0.132691
7 ==> 0.345963% value: 0.0751258
8 ==> 0.0323133% value: 0.00701683
9 ==> 0.0081556% value: 0.00177099
10 ==> 7.44796e-05% value: 1.61732e-05
11 ==> 3.20987e-05% value: 6.97024e-06
12 ==> 2.26398e-05% value: 4.91622e-06
13 ==> 3.67272e-07% value: 7.97531e-08
14 ==> 3.92246e-09% value: 8.51763e-10
15 ==> 3.79906e-09% value: 8.24967e-10
16 ==> 3.63372e-09% value: 7.89062e-10
17 ==> 3.55366e-09% value: 7.71677e-10
18 ==> 3.17073e-09% value: 6.88524e-10
PARAMDIM eigenvalues: 99.0024
5 X 18
Compute PCA constants 
Write constant to file


$ ./fitpca -l
050531-060727-071027-081412-091812-102212 bankstub.txt 
Reading data from bankstub.txt file 
file has 57513 line 
file has 7189 entries 
Read constant from files (c.bin and q.bin)
Using subladder 050531-060727-071027-081412-091812-102212
 numevt: 566
For oneoverpt error 0.397443 0.383462
For phi error 0.0308142 0.0254996
For d0 error 16.2527 4.4808
For cot(tetha/2) error 16.7014 0.978
For z0 error 12.0398 7.27

