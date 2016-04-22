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

> etaval_fit_mean.txt 
> z0val_fit_mean.txt 
> ptval_fit_mean.txt 
> phival_fit_mean.txt 
> etaval_fit_sigma.txt 
> z0val_fit_sigma.txt 
> ptval_fit_sigma.txt 
> phival_fit_sigma.txt 
 

for name in 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
do 
  echo $name >> etaval_fit_mean.txt 
  echo $name >> z0val_fit_mean.txt 
  echo $name >> ptval_fit_mean.txt 
  echo $name >> phival_fit_mean.txt 
  echo $name >> etaval_fit_sigma.txt 
  echo $name >> z0val_fit_sigma.txt 
  echo $name >> ptval_fit_sigma.txt 
  echo $name >> phival_fit_sigma.txt 

  cd $name 
  grep "" output* >> ../etaval.txt 
  grep "" output* >> ../z0val.txt
  grep "q/pt fitted mean" output* >> ../ptval_fit_mean.txt
  grep "q/pt fitted sigma" output* >> ../ptval_fit_sigma.txt
  grep "Phi fitted mean" output* >> ../phival_fit_mean.txt
  grep "Phi fitted sigma" output* >> ../phival_fit_sigma.txt

  grep "Eta fitted mean" output* >> ../etaval_fit_mean.txt
  grep "Eta fitted sigma" output* >> ../etaval_fit_sigma.txt
  grep "z0 fitted mean" output* >> ../z0val_fit_mean.txt
  grep "z0 fitted sigma" output* >> ../z0val_fit_sigma.txt

  cd .. 
done 

