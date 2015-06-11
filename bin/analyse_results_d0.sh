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

tail -n +2 $1 | awk '{print $8}' > d0_orig
tail -n +2 $1 | awk '{print $9}' > d0_pred
python histo.py d0_orig $VALDPT >  d0_orig.hist
python histo.py d0_pred $VALDPT >  d0_pred.hist
tail -n +2 $1 | awk '{print $10}' > d0_diff
python histo.py d0_diff $VALDPT > d0_diff.hist
