/*`
C++ implementation of the PCA fitter

Loriano Storchi: 2016
*/

#include "../interface/PCATrackFitter.h"

#define SF1 1000.0
#define SF2 10.0

namespace 
{
  const int mult_factor = 1e6;
  const int const_mult_factor = mult_factor*1024;
  const int chisq_mult_factor = 1e4;
  const int chisq_const_mult_factor = chisq_mult_factor*1024;
  const int const_w = 25;
  const int add_const_w = 36;


  enum HW_SIGN_TYPE {UNSIGNED, SIGNED};
  /* Function which simulate the HardWare representation of the values : 
   * manage UNSIGNED and SIGNED (2's complement) overflows and accuracy 
   * according to the available dynamic of the binary word from TCB 
   * S.Viret, G.Galbit */
  double binning(double fNumber, int nMSBpowOfTwo, int nBits, HW_SIGN_TYPE signType)
  {
    if (signType == UNSIGNED && fNumber < 0)
      fNumber = -fNumber;
  
    int nLSBpowOfTwo;
	
    //Process the power of two of the LSB for the binary representation
    if (signType == UNSIGNED)
    {
      nLSBpowOfTwo = nMSBpowOfTwo - (nBits-1);
    }	
    else
    {
      //If SIGNED, 1 bit is used for the sign
      nLSBpowOfTwo = nMSBpowOfTwo - (nBits-2);		
    }

    /* Accuracy Simulation */
    //Divide the number by the power of two of the LSB => 
    //the integer part of the new number is the value we are looking for
    fNumber = fNumber / pow(2, nLSBpowOfTwo);
	
    //Remove the fractionnal part by rounding down (for both positive and 
    //negative values), this simulate the HW truncature
    fNumber = floor(fNumber);
	
    //Multiply the number by the power of two of the LSB to get the correct float value
    fNumber = fNumber * pow(2, nLSBpowOfTwo);

    double fBinnedNumber = fNumber;
    /* Overflow Simulation */

    if (signType == UNSIGNED)
    {		  
      //If the number is in UNSIGNED representation			
      fNumber = fmod(fNumber, pow(2, nMSBpowOfTwo+1));
    }
    else
    {		  
      //If the number is in SIGNED representation (2's complement)
      double fTempResult = fNumber - pow(2, nMSBpowOfTwo+1); //substract the possible range to the number

      if (fTempResult >= 0)
      {
        //If there is an overflow, it's a positive one
        fNumber = fmod(fTempResult, pow(2, nMSBpowOfTwo+2)) - pow(2, nMSBpowOfTwo+1);
      }
      else
      {
        //If there is an overflow, it's a negative one (2's complement 
        //format has an asymetric range for positive and negative values)
        fNumber = fmod(fTempResult + pow(2, nLSBpowOfTwo), pow(2, nMSBpowOfTwo+2)) 
          - pow(2, nLSBpowOfTwo) + pow(2, nMSBpowOfTwo+1);
      }  
    }

    //If the new number is different from the previous one, an HW overflow occured
    if (fNumber != fBinnedNumber)
    {
      std::cout << "WARNING HW overflow for the value : " << fBinnedNumber <<
        " resulting value : " << fNumber << " (diff= " << fBinnedNumber-fNumber
        << ")" << std::endl;
    }
	
    return fNumber;
  }

  /*
  #define LAYIDDIM 6
  int stdlayersid[LAYIDDIM] = {5, 6, 7, 8, 9, 10};
  */

  /* TODO to be merged using the  TypeIs... struct */

  /* quick and very dirty */
  template <typename T>
  void dump_element (const pca::matrixpcaconst<T> & in, std::ostream & out)
  {
    out.precision(6);
    for (unsigned int i = 0; i<in.n_rows(); ++i)
    {
      for (unsigned int j = 0; j<in.n_cols(); ++j)
        out << std::scientific << in(i, j) << " ";
      out << std::endl;
    }
  }

  bool hits_to_zrpmatrix_integer (double ci, double si,
      const std::vector<Hit*> & hits, 
      pca::matrixpcaconst<long long int> & zrv, 
      pca::matrixpcaconst<long long int> & phirv, 
      std::string & layersid, 
      std::string & pslayersid,
      int tow)
  {
    std::ostringstream osss, psosss;
    int counter = 0;
    for (unsigned int idx=0; idx<hits.size(); idx++) 
    {
      //double xi = hits[idx]->getX();
      //double yi = hits[idx]->getY();

      double xi = hits[idx]->getX();
      double yi = hits[idx]->getY();

      if ((tow == 19) || (tow == 20) ||
          (tow == 27) || (tow == 28))
      {
        xi = hits[idx]->getX() * ci - hits[idx]->getY() * si;
        yi = hits[idx]->getX() * si + hits[idx]->getY() * ci;
      }

      /* use binning TCB function 
      long long int zi = (long long int) SF1 * hits[idx]->getZ();
      long long int ri = (long long int) SF1 * sqrt(xi*xi+yi*yi);
      long long int pi = SF2 * atan2(yi,xi);
      */
      
      /* as suggested by Geoffrey */
      // Store the integer value of Z (unit = 2^-8 .cm)
      long long int zi =  pow(2, 8) * binning(hits[idx]->getZ(), 8, 18, SIGNED); 
      
      // Store the integer value of R (unit = 2^-10 .cm)
      long long int ri =  pow(2, 10) * binning(sqrt(xi*xi+yi*yi), 6, 18, SIGNED);
      
      // Store the integer value of Phi (unit = 2^-12 .radian)
      long long int pi = pow(2, 16) * binning(atan2(yi,xi), 4, 18, SIGNED);

      zrv(0, counter) = zi;
      phirv(0, counter) = pi;
      ++counter;
      zrv(0, counter) = ri;
      phirv(0, counter) = ri;
      ++counter;

      int lid = (int) hits[idx]->getLayer();

      osss << lid << ":";
      if (lid <= 7)
        psosss << (int) hits[idx]->getLayer() << ":";
    }

    layersid = osss.str();
    layersid.erase(layersid.end()-1);

    pslayersid = psosss.str();
    pslayersid.erase(pslayersid.end()-1);

    return true;
  }

  bool hits_to_zrpmatrix (double ci, double si, 
      const std::vector<Hit*> & hits, 
      pca::matrixpcaconst<double> & zrv, 
      pca::matrixpcaconst<double> & phirv, 
      std::string & layersid, 
      std::string & pslayersid, 
      int tow)
  {
    //std::vector<int> layerids;
    //std::vector<double> riv, piv, ziv;
    
    std::ostringstream osss, psosss;
    int counter = 0;
    for (unsigned int idx=0; idx<hits.size(); ++idx) 
    {
      double xi = hits[idx]->getX();
      double yi = hits[idx]->getY();
 
      // TODO double check this
      if ((tow == 19) || (tow == 20) ||
          (tow == 27) || (tow == 28))
      {
        xi = hits[idx]->getX() * ci - hits[idx]->getY() * si;
        yi = hits[idx]->getX() * si + hits[idx]->getY() * ci;
      }

      //double xi = hits[idx]->getX();
      //double yi = hits[idx]->getY();

      double zi = hits[idx]->getZ();
      double ri = sqrt(xi*xi+yi*yi);
      double pi = atan2(yi,xi);

      /*
      std::cout << "Layer  " << (int) hits[idx]->getLayer() << std::endl;
      std::cout << "PCAxyz " << xi << " " << yi << " " << zi << std::endl;
      std::cout << "PCArpz " << ri << " " << pi << " " << zi << std::endl;

      layerids.push_back((int)hits[idx]->getLayer());
      riv.push_back(ri);
      piv.push_back(pi);
      ziv.push_back(zi);
      */

      zrv(0, counter) = zi;
      phirv(0, counter) = pi;
      ++counter;
      zrv(0, counter) = ri;
      phirv(0, counter) = ri;
      ++counter;

      int lid = (int) hits[idx]->getLayer();

      osss << lid << ":";
      if (lid <= 7)
        psosss << (int) hits[idx]->getLayer() << ":";
    }

    layersid = osss.str();
    layersid.erase(layersid.end()-1);

    pslayersid = psosss.str();
    pslayersid.erase(pslayersid.end()-1);

    //double z1 = 0.0, z3 = 0.0, r1 = 0.0, r3 = 0.0;
    //readoder hits quick test

    /*
    counter = 0;
    for (int i=0; i<LAYIDDIM; ++i)
    {
      std::vector<int>::iterator it = std::find(layerids.begin(),
          layerids.end(),stdlayersid[i]);
      if (it == layerids.end())
        return false;
    
      int pos = (int) std::distance(layerids.begin(), it);
      zrv(0, counter) = ziv[pos];
      phirv(0, counter) = piv[pos];
      ++counter;
      zrv(0, counter) = riv[pos];
      phirv(0, counter) = riv[pos];
      ++counter;

      if (stdlayersid[i] == 5)
      { 
        z1 = ziv[pos];
        r1 = riv[pos];
      }
      else if (stdlayersid[i] == 7)
      {
        z3 = ziv[pos];
        r3 = riv[pos];
      }
    
      std::cout << "Layer  " << layerids[pos] << std::endl;
      std::cout << "PCArpz " << riv[pos] << " " << piv[pos] << " " << ziv[pos] << std::endl;
    }

    double slope = (z3-z1)/(r3-r1);
    double intercept;
    intercept = z1 - slope*r1;

    std::cout << "Intercept est z0: " << intercept << std::endl;
    */
 
    return true;
  }
}

PCATrackFitter::PCATrackFitter():TrackFitter(0)
{
  initialize();
}

PCATrackFitter::PCATrackFitter(int nb):TrackFitter(nb)
{
  initialize();
}

PCATrackFitter::~PCATrackFitter()
{
}

void PCATrackFitter::initialize()
{
  pcacontvct_float_.clear();
  pcacontvct_integer_.clear();
  cleanChi2();
  useinteger_ = false;
  track_ = NULL;
}

void PCATrackFitter::cleanChi2()
{
  chi2vf_.clear();
  chi2vi_.clear();
}

void PCATrackFitter::mergePatterns()
{
}

// Tentative duplicate removal
// Not used for the moment

void PCATrackFitter::mergeTracks()
{
  unsigned int index = 0;
  vector<Track*>::iterator it = tracks.begin();
  
  while(it!=tracks.end())
  {
    Track* newTrack = *it;
    bool found = false;
    for(unsigned int i=0;i<index;i++)
    {
      Track* ref = tracks[i];
      float dpt,dphi,dz,deta;
      dpt = fabs(newTrack->getCurve()-ref->getCurve());
      dphi = fabs(newTrack->getPhi0()-ref->getPhi0());
      dz = fabs(newTrack->getZ0()-ref->getZ0());
      deta = fabs(newTrack->getEta0()-ref->getEta0());
      found = (deta<0.02) &&
	(dphi<0.005) &&
	(dpt<0.1) &&
	(dz<0.3);
      
      if(found)
	break;
    }

    if(found)
      tracks.erase(it);
    else
    {
      index++;
      it++;
    }
  }
}

void PCATrackFitter::setTrack(Track* intc)
{
  track_ = intc;
}

void PCATrackFitter::fit(vector<Hit*> hits)
{
  if (useinteger_)
    this->fit_integer (hits);
  else
    this->fit_float (hits);
}

void PCATrackFitter::fit()
{

  vector<Hit*> activatedHits;

  //////// Get the list of unique stubs from the tracks ///////////

  // TODO not sure will work as it is now need to deal with given 
  //      tracks
  
  set<int> ids;
  int total=0;
  for (unsigned int i=0; i<patterns.size(); ++i)
  {
    vector<Hit*> allHits = patterns[i]->getHits();
    total+=allHits.size();
    for(unsigned int j=0; j<allHits.size(); ++j)
    {
      pair<set<int>::iterator,bool> result = ids.insert(allHits[j]->getID());
      if (result.second == true )
        activatedHits.push_back(allHits[j]);                                
    }
  }
  
  fit(activatedHits);
}

TrackFitter* PCATrackFitter::clone()
{
  PCATrackFitter* fit = new PCATrackFitter(nb_layers);

  // TODO

  return fit;
}

///////////////////////////////////////////////////////////////////////////////
//                                PRIVATE                                    //
///////////////////////////////////////////////////////////////////////////////

void PCATrackFitter::read_float_const_filename (const std::string & in)
{
  std::cout << "Reading " << in << std::endl;
  if (!pca::read_pcaconst_from_file (pcacontvct_float_, in.c_str()))
  {
    std::cerr << "Error while reading constant from " << in << " read only " << 
      pcacontvct_float_.size() << std::endl;
    return;
  }

  std::vector<pca::matrixpcaconst<double> >::const_iterator it = 
      pcacontvct_float_.begin();
  for (; it != pcacontvct_float_.end(); ++it)
  {
    if (it->get_ttype() != pca::FLOATPT)
    {
      std::cerr << "Wrong PCAconst type " << std::endl;
      return;
    }
  }
  std::cout << "Done " << pcacontvct_float_.size() << std::endl;
 
}

void PCATrackFitter::read_integegr_const_filename (const std::string & in)
{
  std::cout << "Reading " << in << std::endl;
  if (!pca::read_pcaconst_from_file (pcacontvct_integer_, in.c_str()))
  {
    std::cerr << "Error while reading constant from " << in << " read only " << 
      pcacontvct_integer_.size() << std::endl;
    return;
  }

  std::vector<pca::matrixpcaconst<long long int> >::const_iterator it = 
      pcacontvct_integer_.begin();
  for (; it != pcacontvct_integer_.end(); ++it)
  {
    if (it->get_ttype() != pca::INTEGPT)
    {
      std::cerr << "Wrong PCAconst type for " << it->get_towerid() << 
        " input file " <<  in << std::endl;
      return;
    }
  }
  std::cout << "Done" << std::endl;
}

void PCATrackFitter::fit_integer(vector<Hit*> hits)
{
  if (pcacontvct_integer_.size() == 0)
  {
    std::cerr << "error PCA const is empty" << std::endl;
    return;
  }

  int tow = sector_id; 

  double sec_phi = (tow%8) * M_PI / 4.0 - 0.4;
  double ci = cos(-sec_phi);
  double si = sin(-sec_phi);

  int charge = +1; /* Try to use only + muons const */ 
  if (track_->getCharge() < 0.0)
    charge = -1;

  std::cout << "Charge from TCB: " << track_->getCharge() << std::endl;

  if (hits.size() == 6)
  {
    pca::matrixpcaconst<long long int> zrv(1, 12), phirv(1, 12);
    std::string layersid, pslayersid;

    if (hits_to_zrpmatrix_integer (ci, si, hits, zrv, phirv, 
          layersid, pslayersid, tow))
    {
      double pt_est = track_->getCurve();
      double eta_est = track_->getEta0();
      double z0_est = track_->getZ0();
      double phi_est = track_->getPhi0();
 
      pca::matrixpcaconst<long long int> cmtx_rz(0, 0);
      pca::matrixpcaconst<long long int> qvec_rz(0, 0); 
      pca::matrixpcaconst<long long int> amtx_rz(0, 0); 
      pca::matrixpcaconst<long long int> kvec_rz(0, 0); 
      pca::matrixpcaconst<long long int> cmtx_rphi(0, 0); 
      pca::matrixpcaconst<long long int> qvec_rphi(0, 0); 
      pca::matrixpcaconst<long long int> amtx_rphi(0, 0); 
      pca::matrixpcaconst<long long int> kvec_rphi(0, 0); 
      
      if (pca::import_pca_const (pcacontvct_integer_, 
                            cmtx_rz, 
                            qvec_rz, 
                            amtx_rz, 
                            kvec_rz, 
                            cmtx_rphi, 
                            qvec_rphi, 
                            amtx_rphi, 
                            kvec_rphi, 
                            eta_est, 
                            pt_est, 
                            charge,
                            layersid, 
                            pslayersid, 
                            tow, 
                            pca::INTEGPT))
      {
        long long int cottheta = 0;
        long long int z0 = 0;
        
        cottheta = qvec_rz(0,0);
        z0 = qvec_rz(0,1);
        for (int i=0; i<(int)cmtx_rz.n_cols(); ++i)
        {
          cottheta += cmtx_rz(0, i) * zrv(0, i);
          z0 += cmtx_rz(1, i) * zrv(0, i);
        }

        double d_cottheta = (double) ((double)cottheta / (double) const_mult_factor);
        double d_z0 = (double) ((double)z0 / (double)const_mult_factor);

        double eta = 0.0e0;
        double theta = atan(1.0e0 / d_cottheta); 
        double tantheta2 = tan (theta/2.0e0); 
        if (tantheta2 < 0.0)
          eta = 1.0e0 * log (-1.0e0 * tantheta2);
        else
          eta = -1.0e0 * log (tantheta2);
        
        long long int coverpt = 0; // pt
        long long int phi = 0;
        
        coverpt = qvec_rphi(0,0);
        phi = qvec_rphi(0,1);
        for (int i=0; i<(int)cmtx_rphi.n_cols(); ++i)
        {
          coverpt += cmtx_rphi(0, i) * phirv(0, i);
          phi += cmtx_rphi(1, i) * phirv(0, i);
        }

        double d_pt = (double) charge / ((double) coverpt / (double) const_mult_factor);
        double d_phi = (double) ((double) phi / (double) const_mult_factor);

        std::cout << " 6oof6 pt:      " << d_pt << " " << pt_est << std::endl;
        std::cout << " 6oof6 phi:     " << d_phi << " " << phi_est << std::endl; 
        std::cout << " 6oof6 eta:     " << eta << " " << eta_est << std::endl;
        std::cout << " 6oof6 z0:      " << d_z0 << " " << z0_est << std::endl;
      }
      else 
      {
        std::cerr << "error while reading PCA const" << std::endl;
      }
    } 
    else
    {
      std::cerr << "error in coord conv" << std::endl;
    }
  }
  else if (hits.size() == 5)
  {
    pca::matrixpcaconst<long long int> zrv(1, 10), phirv(1, 10);
    std::string layersid, pslayersid;

    if (hits_to_zrpmatrix_integer (ci, si, hits, zrv, phirv, 
          layersid, pslayersid, tow))
    {
      double pt_est = track_->getCurve();
      double eta_est = track_->getEta0();
      double z0_est = track_->getZ0();
      double phi_est = track_->getPhi0();
      
      pca::matrixpcaconst<long long int> cmtx_rz(0, 0);
      pca::matrixpcaconst<long long int> qvec_rz(0, 0); 
      pca::matrixpcaconst<long long int> amtx_rz(0, 0); 
      pca::matrixpcaconst<long long int> kvec_rz(0, 0); 
      pca::matrixpcaconst<long long int> cmtx_rphi(0, 0); 
      pca::matrixpcaconst<long long int> qvec_rphi(0, 0); 
      pca::matrixpcaconst<long long int> amtx_rphi(0, 0); 
      pca::matrixpcaconst<long long int> kvec_rphi(0, 0); 
      
      if (import_pca_const (pcacontvct_integer_, 
                            cmtx_rz, 
                            qvec_rz, 
                            amtx_rz, 
                            kvec_rz, 
                            cmtx_rphi, 
                            qvec_rphi, 
                            amtx_rphi, 
                            kvec_rphi, 
                            eta_est, 
                            pt_est, 
                            charge,
                            layersid, 
                            pslayersid, 
                            tow, 
                            pca::INTEGPT))
      {
        long long int cottheta = 0; 
        long long int z0 = 0;
        
        cottheta = qvec_rz(0,0);
        z0 = qvec_rz(0,1);
        for (int i=0; i<(int)cmtx_rz.n_cols(); ++i)
        {
          cottheta += cmtx_rz(0, i) * zrv(0, i);
          z0 += cmtx_rz(1, i) * zrv(0, i);
        }

        double d_cottheta = (double) ((double) cottheta / (double)const_mult_factor);
        double d_z0 = (double) ((double)z0 / (double)const_mult_factor);

        double eta = 0.0e0;
        double theta = atan(1.0e0 / d_cottheta); 
        double tantheta2 = tan (theta/2.0e0); 
        if (tantheta2 < 0.0)
          eta = 1.0e0 * log (-1.0e0 * tantheta2);
        else
          eta = -1.0e0 * log (tantheta2);
        
        long long int coverpt = 0.0; 
        long long int phi = 0.0;
        
        coverpt = qvec_rphi(0,0);
        phi = qvec_rphi(0,1);
        for (int i=0; i<(int)cmtx_rphi.n_cols(); ++i)
        {
          coverpt += cmtx_rphi(0, i) * phirv(0, i);
          phi += cmtx_rphi(1, i) * phirv(0, i);
        }

        double d_pt = (double) charge / ((double) coverpt / (double) const_mult_factor);
        double d_phi = (double) ((double) phi / (double) const_mult_factor);

        std::cout << " 5oof6 pt:      " << d_pt << " " << pt_est << std::endl;
        std::cout << " 5oof6 phi:     " << d_phi << " " << phi_est << std::endl; 
        std::cout << " 5oof6 eta:     " << eta << " " << eta_est << std::endl;
        std::cout << " 5oof6 z0:      " << d_z0 << " " << z0_est << std::endl;
      }
      else 
      {
        std::cerr << "error while reading PCA const" << std::endl;
      }
    } 
    else
    {
      std::cerr << "error in coord conv" << std::endl;
    }

    // TODO 
  }
  else 
  {
    std::cerr << "Stub size ne 5 o 6" << std::endl;
  }

}
 

void PCATrackFitter::fit_float(vector<Hit*> hits)
{
  if (pcacontvct_float_.size() == 0)
  {
    std::cerr << "error PCA const is empty" << std::endl;
    return;
  }

  int tow = sector_id; // The tower ID, necessary to get the phi shift

  //  hits.size() << std::endl;
  
  double sec_phi = (tow%8) * M_PI / 4.0 - 0.4;
  double ci = cos(-sec_phi);
  double si = sin(-sec_phi);

  //std::cout << "tracks_.size() : " << tracks_.size() << std::endl;
  
  /*
  std::ofstream xyzfile;
  xyzfile.open ("PCAxyzfile.txt");

  std::ofstream rpzfile;
  rpzfile.open ("PCArphizfile.txt");
  */

  int charge = +1; /* Try to use only + muons const */ 
  if (track_->getCharge() < 0.0)
    charge = -1;

  std::cout << "Charge from TCB: " << track_->getCharge() << std::endl;

  // Check the charge TODO
  //charge = -1 * charge;

  if (hits.size() == 6)
  {
    pca::matrixpcaconst<double> zrv(1, 12), phirv(1, 12);
    std::string layersid, pslayersid;

    if (hits_to_zrpmatrix (ci, si, hits, zrv, phirv, 
          layersid, pslayersid, tow))
    {
      double pt_est = track_->getCurve();
      double eta_est = track_->getEta0();
      double z0_est = track_->getZ0();
      double phi_est = track_->getPhi0();
      
      /* TEST 
      zrv(0, 0) = -16.1314; zrv(0, 1) = 23.3135;
      zrv(0, 2) = -22.7054; zrv(0, 3) = 34.8367;
      zrv(0, 4) = -32.6362; zrv(0, 5) = 51.9131;
      eta_est = -0.580448;
      z0_est = -2.65203;
      
      phirv(0, 0) = 2.74588; phirv(0, 1) = 23.2105;
      phirv(0, 2) = 2.76857; phirv(0, 3) = 37.0286;
      phirv(0, 4) = 2.79267; phirv(0, 5) = 51.2692;
      phirv(0, 6) = 2.82109; phirv(0, 7) = 67.7289;
      phirv(0, 8) = 2.85643; phirv(0, 9) = 87.8179;
      phirv(0, 10) = 2.89452; phirv(0, 11) = 109.012;
      pt_est = 1.0/0.315806;
      phi_est = 2.70006;
      charge = +1;
      */
      
      pca::matrixpcaconst<double> cmtx_rz(0, 0);
      pca::matrixpcaconst<double> qvec_rz(0, 0); 
      pca::matrixpcaconst<double> amtx_rz(0, 0); 
      pca::matrixpcaconst<double> kvec_rz(0, 0); 
      pca::matrixpcaconst<double> cmtx_rphi(0, 0); 
      pca::matrixpcaconst<double> qvec_rphi(0, 0); 
      pca::matrixpcaconst<double> amtx_rphi(0, 0); 
      pca::matrixpcaconst<double> kvec_rphi(0, 0); 
      
      if (pca::import_pca_const (pcacontvct_float_, 
                            cmtx_rz, 
                            qvec_rz, 
                            amtx_rz, 
                            kvec_rz, 
                            cmtx_rphi, 
                            qvec_rphi, 
                            amtx_rphi, 
                            kvec_rphi, 
                            eta_est, 
                            pt_est, 
                            charge,
                            layersid, 
                            pslayersid, 
                            tow, 
                            pca::FLOATPT))
      {
        /*
        std::cout << "CMTX RZ: " << std::endl;
        dump_element(cmtx_rz, std::cout);
      
        std::cout << "QVEC RZ: " << std::endl;
        dump_element(qvec_rz, std::cout);
      
        std::cout << "CMTX RPHI: " << std::endl;
        dump_element(cmtx_rphi, std::cout);
      
        std::cout << "QVEC RPHI: " << std::endl;
        dump_element(qvec_rphi, std::cout);
        */

        double cottheta = 0.0; // eta
        double z0 = 0.0;
        
        cottheta = qvec_rz(0,0);
        z0 = qvec_rz(0,1);
        for (int i=0; i<(int)cmtx_rz.n_cols(); ++i)
        {
          cottheta += cmtx_rz(0, i) * zrv(0, i);
          z0 += cmtx_rz(1, i) * zrv(0, i);
        }
        
        double coverpt = 0.0; // pt
        double phi = 0.0;
        
        coverpt = qvec_rphi(0,0);
        phi = qvec_rphi(0,1);
        for (int i=0; i<(int)cmtx_rphi.n_cols(); ++i)
        {
          coverpt += cmtx_rphi(0, i) * phirv(0, i);
          phi += cmtx_rphi(1, i) * phirv(0, i);
        }

        if ((tow == 19) || (tow == 20) ||
            (tow == 27) || (tow == 28))
          phi -= sec_phi;

        phi = fmod(phi + M_PI, 2 * M_PI) - M_PI;

        double pt = (double)(charge)/coverpt;
        
        // TODO: checkit theta to eta 
        double eta = 0.0e0;
        double theta = atan(1.0e0 / cottheta); 
        double tantheta2 = tan (theta/2.0e0); 
        if (tantheta2 < 0.0)
          eta = 1.0e0 * log (-1.0e0 * tantheta2);
        else
          eta = -1.0e0 * log (tantheta2);

        int coordim = 6, paramdim = 2;
        double chi2rz = 0.0;
        for (int i=0; i<coordim-paramdim; ++i)
        {
          double val = 0.0;
                                      
          for (int j=0; j<coordim; ++j)
            val += amtx_rz(i,j) * zrv(0, j);

          val -= kvec_rz(0, i);
          
          chi2rz += val*val;
        }

        coordim = 12, paramdim = 2;
        double chi2rphi = 0.0;
        for (int i=0; i<coordim-paramdim; ++i)
        {
          double val = 0.0;
                                      
          for (int j=0; j<coordim; ++j)
            val += amtx_rphi(i,j) * phirv(0, j);

          val -= kvec_rphi(0, i);
          
          chi2rphi += val*val;
        }
        
        std::cout << " 6oof6 pt:      " << pt << " " << pt_est << std::endl;
        std::cout << " 6oof6 phi:     " << phi << " " << phi_est << std::endl; 
        std::cout << " 6oof6 eta:     " << eta << " " << eta_est << std::endl;
        std::cout << " 6oof6 z0:      " << z0 << " " << z0_est << std::endl;
        std::cout << " 6oof6 chirz:   " << chi2rz/4.0 << std::endl;
        std::cout << " 6oof6 chirphi: " << chi2rphi/10.0 << std::endl;

        Track* fit_track = new Track();
        
        fit_track->setCurve(pt);
        fit_track->setPhi0(phi);
        fit_track->setEta0(eta);
        fit_track->setZ0(z0);
                        
        for(unsigned int idx = 0; idx < hits.size(); ++idx)
          fit_track->addStubIndex(hits[idx]->getID());
        
        // TODO: check NDF (14)
        chi2vf_.push_back((chi2rz+chi2rphi)/14.0);
        tracks.push_back(fit_track);
      }
      else 
      {
        std::cerr << "error while reading PCA const" << std::endl;
      }
    } 
    else
    {
      std::cerr << "error in coord conv" << std::endl;
    }

  }
  else if (hits.size() == 5)
  {
    pca::matrixpcaconst<double> zrv(1, 10), phirv(1, 10);
    std::string layersid, pslayersid;

    if (hits_to_zrpmatrix (ci, si, hits, zrv, phirv, 
          layersid, pslayersid, tow))
    {
      double pt_est = track_->getCurve();
      double eta_est = track_->getEta0();
      double z0_est = track_->getZ0();
      double phi_est = track_->getPhi0();
      
      pca::matrixpcaconst<double> cmtx_rz(0, 0);
      pca::matrixpcaconst<double> qvec_rz(0, 0); 
      pca::matrixpcaconst<double> amtx_rz(0, 0); 
      pca::matrixpcaconst<double> kvec_rz(0, 0); 
      pca::matrixpcaconst<double> cmtx_rphi(0, 0); 
      pca::matrixpcaconst<double> qvec_rphi(0, 0); 
      pca::matrixpcaconst<double> amtx_rphi(0, 0); 
      pca::matrixpcaconst<double> kvec_rphi(0, 0); 
      
      if (import_pca_const (pcacontvct_float_, 
                            cmtx_rz, 
                            qvec_rz, 
                            amtx_rz, 
                            kvec_rz, 
                            cmtx_rphi, 
                            qvec_rphi, 
                            amtx_rphi, 
                            kvec_rphi, 
                            eta_est, 
                            pt_est, 
                            charge,
                            layersid, 
                            pslayersid, 
                            tow, 
                            pca::FLOATPT))
      {
        double cottheta = 0.0; // eta
        double z0 = 0.0;
        
        cottheta = qvec_rz(0,0);
        z0 = qvec_rz(0,1);
        for (int i=0; i<(int)cmtx_rz.n_cols(); ++i)
        {
          cottheta += cmtx_rz(0, i) * zrv(0, i);
          z0 += cmtx_rz(1, i) * zrv(0, i);
        }
        
        double coverpt = 0.0; // pt
        double phi = 0.0;
        
        coverpt = qvec_rphi(0,0);
        phi = qvec_rphi(0,1);
        for (int i=0; i<(int)cmtx_rphi.n_cols(); ++i)
        {
          coverpt += cmtx_rphi(0, i) * phirv(0, i);
          phi += cmtx_rphi(1, i) * phirv(0, i);
        }

        if ((tow == 19) || (tow == 20) ||
            (tow == 27) || (tow == 28))
          phi -= sec_phi;

        phi = fmod(phi + M_PI, 2 * M_PI) - M_PI;
        
        double pt = (double)(charge)/coverpt;
        
        // TODO: checkit theta to eta 
        double eta = 0.0e0;
        double theta = atan(1.0e0 / cottheta); 
        double tantheta2 = tan (theta/2.0e0); 
        if (tantheta2 < 0.0)
          eta = 1.0e0 * log (-1.0e0 * tantheta2);
        else
          eta = -1.0e0 * log (tantheta2);

        int coordim = 4, paramdim = 2;
        double chi2rz = 0.0;
        for (int i=0; i<coordim-paramdim; ++i)
        {
          double val = 0.0;
                                      
          for (int j=0; j<coordim; ++j)
            val += amtx_rz(i,j) * zrv(0, j);

          val -= kvec_rz(0, i);
          
          chi2rz += val*val;
        }

        coordim = 10, paramdim = 2;
        double chi2rphi = 0.0;
        for (int i=0; i<coordim-paramdim; ++i)
        {
          double val = 0.0;
                                      
          for (int j=0; j<coordim; ++j)
            val += amtx_rphi(i,j) * phirv(0, j);

          val -= kvec_rphi(0, i);
          
          chi2rphi += val*val;
        }
        
        std::cout << " 5oof6 pt:      " << pt << " " << pt_est << std::endl;
        std::cout << " 5oof6 phi:     " << phi << " " << phi_est << std::endl; 
        std::cout << " 5oof6 eta:     " << eta << " " << eta_est << std::endl;
        std::cout << " 5oof6 z0:      " << z0 << " " << z0_est << std::endl;
        std::cout << " 5oof6 chirz:   " << chi2rz/2.0 << std::endl;
        std::cout << " 5oof6 chirphi: " << chi2rphi/8.0 << std::endl;

        Track* fit_track = new Track();
        
        fit_track->setCurve(pt);
        fit_track->setPhi0(phi);
        fit_track->setEta0(eta);
        fit_track->setZ0(z0);
                        
        for(unsigned int idx = 0; idx < hits.size(); ++idx)
          fit_track->addStubIndex(hits[idx]->getID());
        
        tracks.push_back(fit_track);
        // TODO: check NDF (10)
        //chi2vf_.push_back((chi2rz+chi2rphi)/(10.0));
        // use only rphi 
        chi2vf_.push_back(chi2rphi/8.0);
      }
      else 
      {
        std::cerr << "error while reading PCA const" << std::endl;
      }
    } 
    else
    {
      std::cerr << "error in coord conv" << std::endl;
    }
  }
  else 
  {
    std::cerr << "Stub size ne 5 o 6" << std::endl;
  }

  /*
  xyzfile.close();
  rpzfile.close();

  std::cout << "Close files" << std::endl;
  */
}

