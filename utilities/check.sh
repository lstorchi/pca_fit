for name in 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 ; do echo "==> "$name; cd $name; grep "Chival" output* | wc -l ; cd ..; done
