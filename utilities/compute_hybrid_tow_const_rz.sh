export TOWERID=8
export ROTTFILENAME=MUBANK_pt2To200_tow8.root

export k=-1.7
for i in `seq 1 26`;
do
 export start=$k
 k=$(echo "$k + 0.05" | bc | awk '{printf "%f", $0}')
 echo "Eta range: " $start " "  $k 
 export rage=$start";"$k
 echo $rage

 ./generatepca_split -x -k -R 1 --rz-plane --pt-range="3.0;200.0" --eta-range="$rage" -D $TOWERID $ROTTFILENAME
 ./fitpca_split_single -c "pca_const.txt" -x -k -R 1 --rz-plane --pt-range="3.0;200.0" --eta-range="$rage" -D $TOWERID $ROTTFILENAME

 rm pca_const.txt

 ./generatepca_split -k -R 1 --rz-plane --pt-range="3.0;200.0" --eta-range="$rage" -D $TOWERID $ROTTFILENAME
 ./fitpca_split_single -c "pca_const.txt" -k -R 1 --rz-plane --pt-range="3.0;200.0" --eta-range="$rage" -D $TOWERID $ROTTFILENAME

 rm pca_const.txt
done
