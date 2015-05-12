if [ -z "$1" ]; then 
  echo "$0 results.file.txt"; 
  exit;
fi

VALDPT=100

tail -n +2 $1 | awk '{print 1.0/$1}' > pt_orig
tail -n +2 $1 | awk '{print 1.0/$2}' > pt_pred
python histo.py pt_orig $VALDPT >  pt_orig.hist
python histo.py pt_pred $VALDPT >  pt_pred.hist
tail -n +2 $1 | awk '{print ((1.0/$1 - 1.0/$2)/(1.0/$1))}' > pt_diff
python histo.py pt_diff $VALDPT > pt_diff.hist

tail -n +2 $1 | awk '{print 1.0/$1, ((1.0/$1 - 1.0/$2)/(1.0/$1))}' > pt_diff_vs_pt

tail -n +2 $1 | awk '{print $4}' > phi_orig
tail -n +2 $1 | awk '{print $5}' > phi_pred
python histo.py phi_orig $VALDPT >  phi_orig.hist
python histo.py phi_pred $VALDPT >  phi_pred.hist
tail -n +2 $1 | awk '{print $6}' > phi_diff
python histo.py phi_diff $VALDPT > phi_diff.hist

tail -n +2 $1 | awk '{print $4, $6}' > phi_diff_vs_phi
tail -n +2 $1 | awk '{print 1.0/$1, $6}' > phi_diff_vs_pt

