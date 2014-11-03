./fitpca -V -s 0505-0608-0711-0816-0921-1026 out.1 > outfit 
grep "pt  cmpt" outfit  | awk '{print $3}'  >  ptfit.txt
python ../rootfilereader/histo.py  ptfit.txt > outpt 

