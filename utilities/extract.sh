> chi2val.txt 
> etaval.txt 
> z0val.txt 
> ptval.txt 
> phival.txt 
 
for name in 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
do 
  echo $name >> chi2val.txt 
  echo $name >> etaval.txt 
  echo $name >> z0val.txt 
  echo $name >> ptval.txt 
  echo $name >> phival.txt 
  cd $name 
  grep "Chival" output* >> ../chi2val.txt
  grep "For eta error" output* >> ../etaval.txt 
  grep "For z0 error" output* >> ../z0val.txt
  grep "For q/pt error" output* | grep "%" >> ../ptval.txt
  grep "For phi error" output* >> ../phival.txt
  cd .. 
done 

#grep "Chival" output*
#grep "For eta error" output* 
#grep "For z0 error" output* 
#grep "For q/pt error" output* | grep "%" 
#grep "For phi error" output*
