http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto620

To set up software environment on lxplus:

cmsrel CMSSW_6_2_0_SLHC27
cd CMSSW_6_2_0_SLHC27/src
cmsenv
git-cms-addpkg IOMC/ParticleGuns
git-cms-addpkg IOPool/Input
git-cms-addpkg L1Trigger/TrackFindingAM
git-clone git://github.com/sviret/HL_LHC.git
mv HL_LHC/DataProduction DataProduction
mv HL_LHC/Extractors Extractors
mv HL_LHC/Utils Utils
git checkout CMSSW_6_2_SLHCDEV_X
cp HL_LHC/IOPool_hack/PoolSource.cc IOPool/Input/src/
scramv1 b -j4

cd ../../

git clone git@github.com:sviret/cmssw.git
cd cmssw 
git checkout TrackFindingAM_SV_141215
cd cmssw/L1Trigger/
cp -rfdav TrackFindingAM CMSSW_6_2_0_SLHC27/src/L1Trigger/
cd CMSSW_6_2_0_SLHC27/src/L1Trigger/
scramv1 b -j4


