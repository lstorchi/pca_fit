export TOWERID=8
export ROTTFILENAME=MUBANK_pt2To200_tow8.root
export OUTFILENAME="output_compute_pca_const.out_7bins"

> $OUTFILENAME

declare -a arr=("3.0;7.0" "7.0;12.0" "12.0;18.0" "18.0;25.0" "25.0;50.0" "50.0;100.0" "100.0;200.0")

for i in "${arr[@]}"
do
  echo $i " Gev mu+" >> $OUTFILENAME
  ./generatepca_split -k -R 1 --rphi-plane --charge-sign=+ --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  ./fitpca_split_single -k -N -c "./pca_const.txt" -R 1 --rphi-plane --charge-sign=+ --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  rm results.txt pca_const.txt
  echo ""

done

export OUTFILENAME="output_compute_pca_const.out_17bins"

> $OUTFILENAME

echo "Using " $ROTTFILENAME " for Tower " $TOWERID " out " $OUTFILENAME > $OUTFILENAME

declare -a arr=("3.0;4.0" "4.0;6.0" "6.0;9.0" "9.0;13.0" "13.0;17.0" "17.0;22.0" "22.0;28.0" "28.0;35.0" "35.0;43.0" "43.0;52.0" "52.0;62.0" "62.0;73.0" "73.0;85.0" "85.0;100.0" "100.0;120.0" "120.0;150.0" "150.0;200.0")


for i in "${arr[@]}"
do
  echo $i " Gev mu+" >> $OUTFILENAME
  ./generatepca_split -k -R 1 --rphi-plane --charge-sign=+ --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  ./fitpca_split_single -k -N -c "./pca_const.txt" -R 1 --rphi-plane --charge-sign=+ --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  rm results.txt pca_const.txt
  echo ""
done

export OUTFILENAME="output_compute_pca_const.out_estasplit"

> $OUTFILENAME

declare -a arr=("3.0;7.0" "7.0;12.0" "12.0;18.0" "18.0;25.0" "25.0;50.0" "50.0;100.0" "100.0;200.0")

for i in "${arr[@]}"
do
  echo $i " Gev mu+" >> $OUTFILENAME
  ./generatepca_split -k -R 1 --eta-range="-1.0;-0.4" --rphi-plane --charge-sign=+ --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  ./fitpca_split_single -k -N --eta-range="-1.0;-0.4" -c "./pca_const.txt" -R 1 --rphi-plane --charge-sign=+ --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  rm results.txt pca_const.txt
  echo ""
done

for i in "${arr[@]}"
do
  echo $i " Gev mu+" >> $OUTFILENAME
  ./generatepca_split -k -R 1 --eta-range="-1.7;-1.0" --rphi-plane --charge-sign=+ --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  ./fitpca_split_single -k -N --eta-range="-1.7;-1.0" -c "./pca_const.txt" -R 1 --rphi-plane --charge-sign=+ --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  rm results.txt pca_const.txt
  echo ""
done


