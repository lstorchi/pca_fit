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

tail -n +2 $1 | awk '{print $2}' > x0_orig
tail -n +2 $1 | awk '{print $3}' > x0_pred
python histo.py x0_orig $VALDPT >  x0_orig.hist
python histo.py x0_pred $VALDPT >  x0_pred.hist
tail -n +2 $1 | awk '{print $4}' > x0_diff
python histo.py x0_diff $VALDPT > x0_diff.hist

tail -n +2 $1 | awk '{print $5}' > y0_orig
tail -n +2 $1 | awk '{print $6}' > y0_pred
python histo.py y0_orig $VALDPT >  y0_orig.hist
python histo.py y0_pred $VALDPT >  y0_pred.hist
tail -n +2 $1 | awk '{print $7}' > y0_diff
python histo.py y0_diff $VALDPT > y0_diff.hist
