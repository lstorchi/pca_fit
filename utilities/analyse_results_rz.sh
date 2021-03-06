if [ -z "$1" ]; then 
  echo "$0 results.file.txt"; 
  exit;
fi

VALDPT=80

if [ -z "$2" ]; then 
  VALDPT=80
else
  VALDPT=$2
fi

tail -n +2 $1 | awk '{print $2}' > eta_orig
tail -n +2 $1 | awk '{print $3}' > eta_pred
python histo.py eta_orig $VALDPT >  eta_orig.hist
python histo.py eta_pred $VALDPT >  eta_pred.hist
tail -n +2 $1 | awk '{print $4}' > eta_diff
python histo.py eta_diff $VALDPT > eta_diff.hist

tail -n +2 $1 | awk '{print $5}' > z0_orig
tail -n +2 $1 | awk '{print $6}' > z0_pred
python histo.py z0_orig $VALDPT >  z0_orig.hist
python histo.py z0_pred $VALDPT >  z0_pred.hist
tail -n +2 $1 | awk '{print $7}' > z0_diff
python histo.py z0_diff $VALDPT > z0_diff.hist
