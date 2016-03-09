awk '{print $1-$2}' $1 > pt.txt
awk '{print $4-$5}' $1 > phi.txt
awk '{print $7-$8}' $1 > eta.txt
awk '{print $10-$11}' $1 > d0.txt
awk '{print $13-$14}' $1 > z0.txt

python ../rootfilereader/histo.py  pt.txt > pt.histo
python ../rootfilereader/histo.py  phi.txt > phi.histo
python ../rootfilereader/histo.py  eta.txt > eta.histo
python ../rootfilereader/histo.py  d0.txt > d0.histo
python ../rootfilereader/histo.py  z0.txt > z0.histo



