if [ -z "$1" ]; then 
  echo "$0 results.file.txt"; 
  exit;
fi

tail -n +2 $1 | awk '{print $1}' > onept_orig
tail -n +2 $1 | awk '{print $2}' > onept_pred
python histo.py onept_orig 80 >  onept_orig.hist
python histo.py onept_pred 80 >  onept_pred.hist

tail -n +2 $1 | awk '{print $4}' > phi_orig
tail -n +2 $1 | awk '{print $5}' > phi_pred
python histo.py phi_orig 80 >  phi_orig.hist
python histo.py phi_pred 80 >  phi_pred.hist

tail -n +2 $1 | awk '{print $7}' > eta_orig
tail -n +2 $1 | awk '{print $8}' > eta_pred
python histo.py eta_orig 80 >  eta_orig.hist
python histo.py eta_pred 80 >  eta_pred.hist

tail -n +2 $1 | awk '{print $10}' > d0_orig
tail -n +2 $1 | awk '{print $11}' > d0_pred
python histo.py d0_orig 80 >  d0_orig.hist
python histo.py d0_pred 80 >  d0_pred.hist

tail -n +2 $1 | awk '{print $13}' > z0_orig
tail -n +2 $1 | awk '{print $14}' > z0_pred
python histo.py z0_orig 80 >  z0_orig.hist
python histo.py z0_pred 80 >  z0_pred.hist

