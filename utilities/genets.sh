export k=$1
for i in `seq 1 20`;
do
  export start=$k
  k=$(echo "$k + 0.05" | bc | awk '{printf "%f", $0}')
  echo "$k - 0.025" | bc | awk '{printf "%f\n", $0}'
  export rage=$start";"$k
done
