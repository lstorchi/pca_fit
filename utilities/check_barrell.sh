if [ "$#" -ne 2 ]; then
  echo $0 " rootfile consfile"
  exit
fi

export ROTTFILENAME=$1
export OUTFILENAME="output_check_pca_const.out"
export CONSTFILE=$2

echo "Using " $ROTTFILENAME " and const file " $CONSTFILE " out " $OUTFILENAME 
> $OUTFILENAME

declare -a arr=("3.0;7.0" "7.0;12.0" "12.0;18.0" "18.0;25.0" "25.0;50.0" "50.0;100.0" "100.0;200.0")

echo "Rphi plane" >> $OUTFILENAME 
for i in "${arr[@]}"
do
  echo $i " Gev mu+" >> $OUTFILENAME
  export FFILE="results_rphi_"$i"_p.txt" 
  FILENAME=${FFILE//;/_}
  echo $FILENAME
  ./fitpca_split -c "$CONSTFILE" -k --rphi-plane --charge-sign=+ --pt-range="$i" $ROTTFILENAME >> $OUTFILENAME
  mv results.txt $FILENAME
  echo ""

  echo $i " Gev mu-" >> $OUTFILENAME
  export FFILE="results_rphi_"$i"_n.txt"
  FILENAME=${FFILE//;/_}
  echo $FILENAME
  ./fitpca_split -c "$CONSTFILE" -k --rphi-plane --charge-sign=- --pt-range="$i" $ROTTFILENAME >> $OUTFILENAME
  mv results.txt $FILENAME
  echo ""
done

echo "RZ plane"

export k=$3
for i in `seq 1 20`;
do
  export start=$k
  k=$(echo "$k + 0.05" | bc | awk '{printf "%f", $0}')
  echo "Eta range: " $start " "  $k >> $OUTFILENAME
  export rage=$start";"$k
  ./fitpca_split -c "$CONSTFILE" -x --rz-plane --eta-range="$rage" -k $ROTTFILENAME >> $OUTFILENAME
  mv results.txt results_rz_$i.txt
  echo ""
done


echo "Rphi plane" >> $OUTFILENAME
echo "Five out of six" >> $OUTFILENAME

for i in "${arr[@]}"
do
  for j in `seq 5 10`
  do 
    echo "Removing layer " $j
    echo $i " Gev mu+" >> $OUTFILENAME
    export FFILE="results_rphi_"$i"_p_fk"$j".txt"
    FILENAME=${FFILE//;/_}
    echo $FILENAME
    ./fitpca_split -c "$CONSTFILE" --fk-five-hits=$j -k --rphi-plane --charge-sign=+ --pt-range="$i" $ROTTFILENAME >> $OUTFILENAME
    mv results.txt $FILENAME
    echo ""

    echo $i " Gev mu-" >> $OUTFILENAME
    export FFILE="results_rphi_"$i"_n_fk"$j".txt"
    FILENAME=${FFILE//;/_}
    echo $FILENAME
    ./fitpca_split -c "$CONSTFILE" --fk-five-hits=$j -k --rphi-plane --charge-sign=- --pt-range="$i" $ROTTFILENAME >> $OUTFILENAME
    mv results.txt $FILENAME
    echo ""
  done

done

echo "RZ plane" >> $OUTFILENAME
echo "Five out of six" >> $OUTFILENAME

export k=$3
for i in `seq 1 20`;
do
  export start=$k
  k=$(echo "$k + 0.05" | bc | awk '{printf "%f", $0}')
  echo "Eta range: " $start " "  $k >> $OUTFILENAME
  export rage=$start";"$k

  ./fitpca_split -c "$CONSTFILE" -x --rz-plane --eta-range="$rage" -k --fk-five-hits=5 $ROTTFILENAME >> $OUTFILENAME
  mv results.txt results_fk5_rz_$i.txt
 
  ./fitpca_split -c "$CONSTFILE" -x --rz-plane --eta-range="$rage" -k --fk-five-hits=6 $ROTTFILENAME >> $OUTFILENAME
  mv results.txt results_fk6_rz_$i.txt
 
  ./fitpca_split -c "$CONSTFILE" -x --rz-plane --eta-range="$rage" -k --fk-five-hits=7 $ROTTFILENAME >> $OUTFILENAME
  mv results.txt results_fk7_rz_$i.txt
done

