
To set up software environment on lxplus:

cmsrel CMSSW_6_2_0_SLHC27
cd CMSSW_6_2_0_SLHC27/src
cmsenv
git-cms-addpkg L1Trigger/TrackFindingAM
git-pull https://github.com/sviret/cmssw L1Tracking_PCA_sviret_040316
scramv1 b -j4

now clone the cpsa_fit repo to use the last CMSSW code version

cd  ../../
git clone https://github.com/lstorchi/pca_fit
cd pca_fit/CMSSW/src/L1Trigger/TrackFindingAM/
cp -rfadv * ../../../../../CMSSW_6_2_0_SLHC27/src/L1Trigger/TrackFindingAM/
cd ../../../../../CMSSW_6_2_0_SLHC27/src/L1Trigger/TrackFindingAM/
scramv1 b -j4


