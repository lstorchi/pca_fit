for num in 8 9 10 11 12 13 14 15
do 
  mkdir $num
  cd $num
  ln -s ../fitpca_split_single .
  ln -s ../generatepca_split . 
  ln -s ../compute_hybrid_tow_const.sh .
  export FFILE="../Hybrid/MUBANK_pt2To200_tow"$num".root"
  ./compute_hybrid_tow_const.sh $FFILE $num -1.4 1> out 2> err &
  cd ../ 
done
