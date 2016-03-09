echo "3 7 Gev"
./generatepca_split -k --rphi-plane --charge-sign=+ --pt-range="3.0;7.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
./fitpca_split -k --rphi-plane --charge-sign=+ --pt-range="3.0;7.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
mv results.txt results_3_7_p.txt
echo ""

echo "7 12 Gev"
./generatepca_split -k --rphi-plane --charge-sign=+ --pt-range="7.0;12.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
./fitpca_split -k --rphi-plane --charge-sign=+ --pt-range="7.0;12.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
mv results.txt results_7_12_p.txt
echo ""

echo "12 18 Gev"
./generatepca_split -k --rphi-plane --charge-sign=+ --pt-range="12.0;18.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
./fitpca_split -k --rphi-plane --charge-sign=+ --pt-range="12.0;18.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
mv results.txt results_12_18_p.txt
echo ""

echo "18 25 Gev"
./generatepca_split -k --rphi-plane --charge-sign=+ --pt-range="18.0;25.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
./fitpca_split -k --rphi-plane --charge-sign=+ --pt-range="18.0;25.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
mv results.txt results_18_25_p.txt
echo ""

echo "25 50 Gev"
./generatepca_split -k --rphi-plane --charge-sign=+ --pt-range="25.0;50.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
./fitpca_split -k --rphi-plane --charge-sign=+ --pt-range="25.0;50.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
mv results.txt results_25_50_p.txt
echo ""

echo "50 100 Gev"
./generatepca_split -k --rphi-plane --charge-sign=+ --pt-range="50.0;100.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
./fitpca_split -k --rphi-plane --charge-sign=+ --pt-range="50.0;100.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
mv results.txt results_50_100_p.txt
echo ""

echo "100 200 Gev"
./generatepca_split -k --rphi-plane --charge-sign=+ --pt-range="100.0;200.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
./fitpca_split -k --rphi-plane --charge-sign=+ --pt-range="100.0;200.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
mv results.txt results_100_200_p.txt
echo ""

echo "3 7 Gev"
./generatepca_split -k --rphi-plane --charge-sign=- --pt-range="3.0;7.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
./fitpca_split -k --rphi-plane --charge-sign=- --pt-range="3.0;7.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
mv results.txt results_3_7_n.txt
echo ""

echo "7 12 Gev"
./generatepca_split -k --rphi-plane --charge-sign=- --pt-range="7.0;12.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
./fitpca_split -k --rphi-plane --charge-sign=- --pt-range="7.0;12.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
mv results.txt results_7_12_n.txt
echo ""

echo "12 18 Gev"
./generatepca_split -k --rphi-plane --charge-sign=- --pt-range="12.0;18.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
./fitpca_split -k --rphi-plane --charge-sign=- --pt-range="12.0;18.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
mv results.txt results_12_18_n.txt
echo ""

echo "18 25 Gev"
./generatepca_split -k --rphi-plane --charge-sign=- --pt-range="18.0;25.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
./fitpca_split -k --rphi-plane --charge-sign=- --pt-range="18.0;25.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
mv results.txt results_18_25_n.txt
echo ""

echo "25 50 Gev"
./generatepca_split -k --rphi-plane --charge-sign=- --pt-range="25.0;50.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
./fitpca_split -k --rphi-plane --charge-sign=- --pt-range="25.0;50.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
mv results.txt results_25_50_n.txt
echo ""

echo "50 100 Gev"
./generatepca_split -k --rphi-plane --charge-sign=- --pt-range="50.0;100.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
./fitpca_split -k --rphi-plane --charge-sign=- --pt-range="50.0;100.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
mv results.txt results_50_100_n.txt
echo ""

echo "100 200 Gev"
./generatepca_split -k --rphi-plane --charge-sign=- --pt-range="100.0;200.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
./fitpca_split -k --rphi-plane --charge-sign=- --pt-range="100.0;200.0" ./MUBANK_pt2To200_phi11To29_etaM06To04.root
mv results.txt results_100_200_n.txt
echo ""

