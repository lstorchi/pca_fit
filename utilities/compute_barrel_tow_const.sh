if [ "$#" -ne 2 ]; then
  echo $0 " rootfile towerid "
  exit
fi

export TOWERID=$2
export ROTTFILENAME=$1
export OUTFILENAME="output_compute_pca_const.out"

echo "Using " $ROTTFILENAME " for Tower " $TOWERID " out " $OUTFILENAME 
> $OUTFILENAME

declare -a arr=("3.0;7.0" "7.0;12.0" "12.0;18.0" "18.0;25.0" "25.0;50.0" "50.0;100.0" "100.0;200.0")

for i in "${arr[@]}"
do
  echo $i " Gev mu+" >> $OUTFILENAME
  export FFILE="results_rphi_"$i"_p.txt" 
  FILENAME=${FFILE//;/_}
  echo $FILENAME
  ./generatepca_split -k --rphi-plane --charge-sign=+ --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  ./fitpca_split -k --rphi-plane --charge-sign=+ --pt-range="3.0;7.0" $ROTTFILENAME >> $OUTFILENAME
  mv results.txt $FILENAME
  echo ""

  echo $i " Gev mu-" >> $OUTFILENAME
  export FFILE="results_rphi_"$i"_n.txt"
  FILENAME=${FFILE//;/_}
  echo $FILENAME
  ./generatepca_split -k --rphi-plane --charge-sign=- --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  ./fitpca_split -k --rphi-plane --charge-sign=- --pt-range="$i" $ROTTFILENAME >> $OUTFILENAME
  mv results.txt $FILENAME
  echo ""

  for j in seq 5 10
  do 
    echo $i " Gev mu+" >> $OUTFILENAME
    export FFILE="results_rphi_"$i"_p_fk"$j".txt"
    FILENAME=${FFILE//;/_}
    echo $FILENAME
    ./generatepca_split --fk-five-hits=$j -k --rphi-plane --charge-sign=+ --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
    ./fitpca_split --fk-five-hits=$j -k --rphi-plane --charge-sign=+ --pt-range="$i" $ROTTFILENAME >> $OUTFILENAME
    mv results.txt $FILENAME
    echo ""

    echo $i " Gev mu-" >> $OUTFILENAME
    export FFILE="results_rphi_"$i"_n_fk"$j".txt"
    FILENAME=${FFILE//;/_}
    echo $FILENAME
    ./generatepca_split --fk-five-hits=$j -k --rphi-plane --charge-sign=- --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
    ./fitpca_split --fk-five-hits=$j -k --rphi-plane --charge-sign=- --pt-range="$i" $ROTTFILENAME >> $OUTFILENAME
    mv results.txt $FILENAME
    echo ""
  done

done

export k="-0.6" 
for i in `seq 1 20`;
do
  export start=$k
  k=$(echo "$k + 0.05" | bc | awk '{printf "%f", $0}')
  echo "Eta range: " $start " "  $k >> $OUTFILENAME
  export rage=$start";"$k
  ./generatepca_split -x --rz-plane --eta-range="$rage" -k --fk-five-hits=5 -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  ./fitpca_split -x --rz-plane --eta-range="$rage" -k --fk-five-hits=5 $ROTTFILENAME >> $OUTFILENAME
  mv results.txt results_fk5_rz_$i.txt
 
  ./generatepca_split -x --rz-plane --eta-range="$rage" -k --fk-five-hits=6 -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  ./fitpca_split -x --rz-plane --eta-range="$rage" -k --fk-five-hits=6 $ROTTFILENAME >> $OUTFILENAME
  mv results.txt results_fk6_rz_$i.txt
 
  ./generatepca_split -x --rz-plane --eta-range="$rage" -k --fk-five-hits=7 -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  ./fitpca_split -x --rz-plane --eta-range="$rage" -k --fk-five-hits=7 $ROTTFILENAME >> $OUTFILENAME
  mv results.txt results_fk7_rz_$i.txt
 
  ./generatepca_split -x --rz-plane --eta-range="$rage" -k -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  ./fitpca_split -x --rz-plane --eta-range="$rage" -k $ROTTFILENAME >> $OUTFILENAME
  mv results.txt results_rz_$i.txt
  echo ""
done

