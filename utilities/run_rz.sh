export k="-0.6" 
for i in `seq 1 20`;
do
  export start=$k
  k=$(echo "$k + 0.05" | bc | awk '{printf "%f", $0}')
  echo "Eta range: " $start " "  $k
  export rage=$start";"$k
  ./generatepca_split -x --rz-plane --eta-range="$rage" -k --fk-five-hits=5 ./MUBANK_pt2To200_phi11To29_etaM06To04.root
  ./fitpca_split -x --rz-plane --eta-range="$rage" -k --fk-five-hits=5  ./MUBANK_pt2To200_phi11To29_etaM06To04.root
  ./generatepca_split -x --rz-plane --eta-range="$rage" -k --fk-five-hits=6 ./MUBANK_pt2To200_phi11To29_etaM06To04.root
  ./fitpca_split -x --rz-plane --eta-range="$rage" -k --fk-five-hits=6  ./MUBANK_pt2To200_phi11To29_etaM06To04.root
  ./generatepca_split -x --rz-plane --eta-range="$rage" -k --fk-five-hits=7 ./MUBANK_pt2To200_phi11To29_etaM06To04.root
  ./fitpca_split -x --rz-plane --eta-range="$rage" -k --fk-five-hits=7  ./MUBANK_pt2To200_phi11To29_etaM06To04.root
  ./generatepca_split -x --rz-plane --eta-range="$rage" -k ./MUBANK_pt2To200_phi11To29_etaM06To04.root
  ./fitpca_split -x --rz-plane --eta-range="$rage" -k ./MUBANK_pt2To200_phi11To29_etaM06To04.root
  mv results.txt results_$i.txt
  echo ""
done

