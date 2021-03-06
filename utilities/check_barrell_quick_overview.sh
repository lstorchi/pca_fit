if [ "$#" -ne 4 ]; then
  echo $0 " rootfile consfile etamin towerid"
  exit
fi

export ROTTFILENAME=$1
export OUTFILENAME="output_check_pca_const.out"
export CONSTFILE=$2
export TOWERID=$4

echo "Using " $ROTTFILENAME " and out " $OUTFILENAME " towerid " $TOWERID
> $OUTFILENAME

#declare -a arr=("3.0;7.0" "7.0;12.0" "12.0;18.0" "18.0;25.0" "25.0;50.0" "50.0;100.0" "100.0;200.0")
declare -a arr=("3.0;7.0")

echo "Rphi plane" >> $OUTFILENAME 
for i in "${arr[@]}"
do
  echo $i " Gev mu+" >> $OUTFILENAME
  export FFILE="results_rphi_"$i"_p.txt" 
  FILENAME=${FFILE//;/_}
  echo $FILENAME
  ./generatepca_split -k --rphi-plane --charge-sign=+ --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  ./fitpca_split -N -k --rphi-plane --charge-sign=+ --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  rm -f results.txt
  echo ""

  echo $i " Gev mu-" >> $OUTFILENAME
  export FFILE="results_rphi_"$i"_n.txt"
  FILENAME=${FFILE//;/_}
  echo $FILENAME
  ./generatepca_split -k --rphi-plane --charge-sign=- --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  ./fitpca_split -N -k --rphi-plane --charge-sign=- --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  rm -f results.txt 
  echo ""
done

echo "RZ plane"

export k=$3
#for i in `seq 1 20`;
for i in `seq 1 1`;
do
  export start=$k
  k=$(echo "$k + 0.05" | bc | awk '{printf "%f", $0}')
  echo "Eta range: " $start " "  $k >> $OUTFILENAME
  export rage=$start";"$k
  ./generatepca_split -x --rz-plane --eta-range="$rage" -k -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  ./fitpca_split -N -x --rz-plane --eta-range="$rage" -k -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  rm -f results.txt
  echo ""
done

