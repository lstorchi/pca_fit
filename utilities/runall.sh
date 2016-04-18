for num in 16 17 18 19 20 21 22 23
do 
  mkdir $num
  cd $num
  ln -s ../fitpca_split .
  ln -s ../generatepca_split . 
  ln -s ../check_barrell_quick_overview.sh .
  export FFILE="/pub/redo/MUBANK/MUBANK_pt2To200_tow"$num".root"
  export CFILE="../../files/barrel_tow"$num"_pca_const.txt"
  ./check_barrell_quick_overview.sh $FFILE $CFILE -0.6 $num 1> out 2> err &
  cd ../ 
done

for num in 24 25 26 27 28 29 30 31
do 
  mkdir $num
  cd $num
  ln -s ../fitpca_split .
  ln -s ../generatepca_split .
  ln -s ../check_barrell_quick_overview.sh .
  export FFILE="/pub/redo/MUBANK/MUBANK_pt2To200_tow"$num".root"
  export CFILE="../../files/barrel_tow"$num"_pca_const.txt"
  ./check_barrell_quick_overview.sh $FFILE $CFILE -0.4 $num 1> out 2> err &
  cd ../ 
done
