declare -a arr=("3.0 7.0" "7.0 12.0" "12.0 18.0" "18.0 25.0" "25.0 50.0" "50.0 100.0" "100.0 200.0")

for i in "${arr[@]}"
do
  echo $i " mu+" 
  echo $i " mu-"
done

