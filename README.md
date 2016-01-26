README (Loriano Storchi) loriano@storchi.org
---------------------------------------------------------------------
- basic pca fitting methods. Here  you  find  the  rootfile  reader 
  as well as the generatepca constants code and the code performing 
  the fitting.
- IMPORTANTE NOTE: This is just a test code for internal use only, 
  nothing more than a basic rapid quick and dirty prototyping (i.e.
  the code should be totally restructured in a proper OO fashion)
- also HEAD rev should be safe enough to be used. Step by step:
- reminder : rel_0_4_0 need to be fixed but is the last tag before 
  starting revision, removing several options (as the one to use a 
  single parameter at time, the one to fit also d0 or x0 in both planes) 
  but need to be fixed start the fixing in the branch

  a) dowload the code  

  $ git clone https://github.com/lstorchi/pca_fit
  $ cd pca_fit
  
  HEAD should be now safe enough to be used

  b) Download the armadillo library from http://arma.sourceforge.net/download.html 
  
  c) Modify config.mk to indicate the correct location of this armadillo library, if needed.
  
  d) Set the environment variable LD_LIBRARY_PATH to point to the armadillo library, id needed.

All codes are provided with some basic help, just use the  --help option:

./bin/readrootfile  read root file

split r-phi and r-z plane:

./bin/generatepca_split generate pca constants.

./bin/fitpca_split perform spitting.

old approach to use specific sub-tower and / or sub-sub-tower has been removed.

Example: 

TODO


