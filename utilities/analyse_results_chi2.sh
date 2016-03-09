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

tail -n +2 $1 | awk '{print $8}' > chi2_10
tail -n +2 $1 | awk '{print $9}' > chi2_11
python histo.py chi2_10 $VALDPT >  chi2_10.hist
python histo.py chi2_11 $VALDPT >  chi2_11.hist
