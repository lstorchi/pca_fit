if [ -z "$1" ]; then 
  echo "$0 results.file.txt"; 
  exit;
fi

VALDPT=80

tail -n +2 $1 | awk '{print $2}' > onept_orig
tail -n +2 $1 | awk '{print $3}' > onept_pred
python histo.py onept_orig $VALDPT >  onept_orig.hist
python histo.py onept_pred $VALDPT >  onept_pred.hist
tail -n +2 $1 | awk '{print $3-$2}' > onept_diff
python histo.py onept_diff $VALDPT > onept_diff.hist

tail -n +2 $1 | awk '{print $5}' > phi_orig
tail -n +2 $1 | awk '{print $6}' > phi_pred
python histo.py phi_orig $VALDPT >  phi_orig.hist
python histo.py phi_pred $VALDPT >  phi_pred.hist
tail -n +2 $1 | awk '{print $6-$5}' > phi_diff
python histo.py phi_diff $VALDPT > phi_diff.hist

tail -n +2 $1 | awk '{print $8}' > eta_orig
tail -n +2 $1 | awk '{print $9}' > eta_pred
python histo.py eta_orig $VALDPT >  eta_orig.hist
python histo.py eta_pred $VALDPT >  eta_pred.hist
tail -n +2 $1 | awk '{print $9-$8}' > eta_diff
python histo.py eta_diff $VALDPT > eta_diff.hist

tail -n +2 $1 | awk '{print $11}' > d0_orig
tail -n +2 $1 | awk '{print $12}' > d0_pred
python histo.py d0_orig $VALDPT >  d0_orig.hist
python histo.py d0_pred $VALDPT >  d0_pred.hist
tail -n +2 $1 | awk '{print $12-$11}' > d0_diff
python histo.py d0_diff $VALDPT > d0_diff.hist

tail -n +2 $1 | awk '{print $14}' > z0_orig
tail -n +2 $1 | awk '{print $15}' > z0_pred
python histo.py z0_orig $VALDPT >  z0_orig.hist
python histo.py z0_pred $VALDPT >  z0_pred.hist
tail -n +2 $1 | awk '{print $15-$14}' > z0_diff
python histo.py z0_diff $VALDPT > z0_diff.hist

