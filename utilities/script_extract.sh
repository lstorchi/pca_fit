  grep "oof6 pt:" AMFIT_test.out  | awk '{print $3}'  > pt_float
  grep "oof6 int pt:" AMFIT_test.out  | awk '{print $4}'  > pt_int
  grep "oof6 phi:" AMFIT_test.out  | awk '{print $3}'  > phi_float
  grep "oof6 int phi:" AMFIT_test.out  | awk '{print $4}'  > phi_int
  grep "oof6 eta:" AMFIT_test.out  | awk '{print $3}'  > est_float
  grep "oof6 int eta:" AMFIT_test.out  | awk '{print $4}'  > eta_int
  grep "oof6 eta:" AMFIT_test.out  | awk '{print $3}'  > eta_float
  grep "oof6 z0:" AMFIT_test.out  | awk '{print $3}'  > z0_float
  grep "oof6 int z0:" AMFIT_test.out  | awk '{print $4}'  > z0_int
  grep "oof6 int chi2rz:" AMFIT_test.out  | awk '{print $4}'  > chi2rz_int
  grep "oof6 chi2rz:" AMFIT_test.out  | awk '{print $3}'  > chi2rz_float
  grep "oof6 chi2rphi:" AMFIT_test.out  | awk '{print $3}'  > chi2rphi_float
  grep "oof6 int chi2rphi:" AMFIT_test.out  | awk '{print $4}'  > chi2rphi_int
  paste phi_float phi_int > phi
  paste pt_float pt_int > pt
  paste eta_float eta_int > eta
  paste z0_float z0_int > z0
  paste chi2rz_float chi2rz_int > chi2rz
  paste chi2rphi_float chi2rphi_int > chi2rphi
  paste chi2rphi_float chi2rphi_int > chi2rphi
  echo "pt"
  python histo_bin.py pt
  echo "phi"
  python histo_bin.py phi
  echo "eta"
  python histo_bin.py eta
  echo "z0"
  python histo_bin.py z0
  echo "chi2rz"
  python histo_bin.py chi2rz
  echo "chi2rphi"
  python histo_bin.py chi2rphi
