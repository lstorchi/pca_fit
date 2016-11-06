if [ "$#" -ne 3 ]; then
  echo $0 " rootfile towerid etaminval"
  exit
fi

export TOWERID=$2
export ROTTFILENAME=$1
export OUTFILENAME="output_compute_pca_const.out"

for seq in $(cat TCB_layer_seq_sec"$TOWERID".txt)
do 
  n=$(grep -o ";" <<< "$seq" | wc -l)
  if [ $n == 5 ]; then
    echo $n " " $seq 
  fi
done

for seq in $(cat TCB_layer_seq_sec"$TOWERID".txt)
do 
  n=$(grep -o ";" <<< "$seq" | wc -l)
  if [ $n == 2 ]; then
    echo $n " " $seq 
  fi
done

echo "Using " $ROTTFILENAME " for Tower " $TOWERID " out " $OUTFILENAME " etamin " $3
> $OUTFILENAME

declare -a arr=("2.0;3.0" "3.0;7.0" "7.0;12.0" "12.0;18.0" "18.0;25.0" "25.0;50.0" "50.0;100.0" "100.0;200.0")

export EXTRAOPT="-X 1000000"
export EXTRAOPT=""

for i in "${arr[@]}"
do
  echo $i " Gev mu+" >> $OUTFILENAME
  ./generatepca_split $EXTRAOPT -R 1 -k --rphi-plane --charge-sign=+ --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  echo ""

  echo $i " Gev mu-" >> $OUTFILENAME
  ./generatepca_split $EXTRAOPT -R 1 -k --rphi-plane --charge-sign=- --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  echo ""
done

export k=$3
for i in `seq 1 26`;
do
  export start=$k
  k=$(echo "$k + 0.05" | bc | awk '{printf "%f", $0}')
  echo "Eta range: " $start " "  $k >> $OUTFILENAME
  export rage=$start";"$k
  ./generatepca_split $EXTRAOPT -R 1 -x --rz-plane --eta-range="$rage" -k -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
  echo ""
done

for i in "${arr[@]}"
do
  for seq in $(cat TCB_layer_seq_sec"$TOWERID".txt)
  do 
    n=$(grep -o ";" <<< "$seq" | wc -l)
    if [ $n == 5 ]; then
      echo $i " Gev mu+" >> $OUTFILENAME
      echo $seq  >> $OUTFILENAME
      ./generatepca_split $EXTRAOPT -R 1 --specific-fk-53-hits="$seq" -k --rphi-plane --charge-sign=+ --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
      echo ""

      echo $i " Gev mu-" >> $OUTFILENAME
      echo $seq  >> $OUTFILENAME
      ./generatepca_split $EXTRAOPT -R 1 --specific-fk-53-hits="$seq" -k --rphi-plane --charge-sign=- --pt-range="$i" -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
      echo ""
    fi
  done
done

export k=$3
for i in `seq 1 26`;
do
  export start=$k
  k=$(echo "$k + 0.05" | bc | awk '{printf "%f", $0}')
  export rage=$start";"$k

  for seq in $(cat TCB_layer_seq_sec"$TOWERID".txt)
  do 
    n=$(grep -o ";" <<< "$seq" | wc -l)
    if [ $n == 2 ]; then
      echo "Eta range: " $start " "  $k >> $OUTFILENAME
      ./generatepca_split $EXTRAOPT -R 1 --specific-fk-53-hits="$seq"-x --rz-plane --eta-range="$rage" -k -D $TOWERID $ROTTFILENAME >> $OUTFILENAME
      echo ""
    fi
  done
done
