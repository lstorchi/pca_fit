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

R-Phi (2.0 <= pt <= 7.0): 

./generatepca_split --rphi-plane --charge-sign=+ --pt-range="2.0;7.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
./fitpca_split --rphi-plane --charge-sign=+ --pt-range="2.0;7.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root

R-Z (D_eta = 0.0.5):

export k="-0.6" 
for i in `seq 1 20`;
do
  export start=$k
  k=$(echo "$k + 0.05" | bc | awk '{printf "%f", $0}')
  echo "Eta range: " $start " "  $k
  export rage=$start";"$k
  ./generatepca_split -x --rz-plane --eta-range="$rage" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
  ./fitpca_split -x --rz-plane --eta-range="$rage" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
  mv results.txt results_$i.txt
  echo ""
done

