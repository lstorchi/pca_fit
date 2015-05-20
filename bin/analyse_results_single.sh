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

tail -n +2 $1 | awk '{print $1}' > single_orig
tail -n +2 $1 | awk '{print $2}' > single_pred
python histo.py single_orig $VALDPT >  single_orig.hist
python histo.py single_pred $VALDPT >  single_pred.hist
tail -n +2 $1 | awk '{print $3}' > single_diff
python histo.py single_diff $VALDPT > single_diff.hist
