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

tail -n +2 $1 | awk '{print $1}' > pt_orig
tail -n +2 $1 | awk '{print $2}' > pt_pred
python histo.py pt_orig $VALDPT >  pt_orig.hist
python histo.py pt_pred $VALDPT >  pt_pred.hist
tail -n +2 $1 | awk '{print $3}' > pt_diff
python histo.py pt_diff $VALDPT > pt_diff.hist

tail -n +2 $1 | awk '{print $4}' > phi_orig
tail -n +2 $1 | awk '{print $5}' > phi_pred
python histo.py phi_orig $VALDPT >  phi_orig.hist
python histo.py phi_pred $VALDPT >  phi_pred.hist
tail -n +2 $1 | awk '{print $6}' > phi_diff
python histo.py phi_diff $VALDPT > phi_diff.hist

