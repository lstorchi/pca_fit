/*
C++ implementation of the PCA fitter

Loriano Storchi: 2016
*/

#include "../interface/PCATrackFitter.h"

namespace 
{
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
      pca::matrixpcaconst<int32_t> & zrv, 
      pca::matrixpcaconst<int32_t> & phirv, 
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

      int32_t zi = (int32_t) hits[idx]->getZ();
      int32_t ri = (int32_t) sqrt(xi*xi+yi*yi);
      int32_t pi = 50000 * atan2(yi,xi);

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

  bool import_pca_const (
      std::vector<pca::matrixpcaconst<int32_t> > & vct,
      pca::matrixpcaconst<int32_t> & cmtx_rz, 
      pca::matrixpcaconst<int32_t> & qvec_rz, 
      pca::matrixpcaconst<int32_t> & amtx_rz, 
      pca::matrixpcaconst<int32_t> & kvec_rz, 
      pca::matrixpcaconst<int32_t> & cmtx_rphi, 
      pca::matrixpcaconst<int32_t> & qvec_rphi, 
      pca::matrixpcaconst<int32_t> & amtx_rphi, 
      pca::matrixpcaconst<int32_t> & kvec_rphi, 
      double eta, double pt, 
      int chargesignin,
      const std::string & layersid,
      const std::string & pslayersid,
      int towerid)
  {
    int hwmanygot = 0, hwmanygotrphi = 0, hwmanygotrz = 0;;
    std::vector<pca::matrixpcaconst<int32_t> >::const_iterator it = 
      vct.begin();
    for (; it != vct.end(); ++it)
    {
      double ptmin, ptmax, etamin, etamax;
      std::string actuallayids;
      int chargesign;
  
      it->get_ptrange(ptmin, ptmax);
      it->get_etarange(etamin, etamax);
      chargesign = it->get_chargesign();
      actuallayids = it->get_layersids();

      if (towerid == it->get_towerid())
      {
        if (it->get_ttype() == pca::matrixpcaconst<int32_t>::INTEGPT)
        {
          if (it->get_plane_type() == pca::matrixpcaconst<int32_t>::RZ)
          {
            if (actuallayids == pslayersid)
            {
              if ((eta >= etamin) && (eta <= etamax)) 
              {
                switch(it->get_const_type())
                {
                  case pca::matrixpcaconst<int32_t>::QVEC :
                    qvec_rz = *it;
                    hwmanygot++;
                    hwmanygotrz++;
                    break;
                  case pca::matrixpcaconst<int32_t>::KVEC :
                    kvec_rz = *it;
                    hwmanygot++;
                    hwmanygotrz++;
                    break;
                  case pca::matrixpcaconst<int32_t>::CMTX :
                    cmtx_rz = *it;
                    hwmanygot++;
                    hwmanygotrz++;
                    break;
                  case pca::matrixpcaconst<int32_t>::AMTX :
                    amtx_rz = *it;
                    hwmanygot++;
                    hwmanygotrz++;
                    break;
                  default:
                    break;
                }
              } 
            }
          }
          else if (it->get_plane_type() == pca::matrixpcaconst<int32_t>::RPHI)
          {
            if (actuallayids == layersid)
            {
              if (chargesignin == chargesign)
              {
                if ((pt >= ptmin) && (pt <= ptmax))
                {
                  switch(it->get_const_type())
                  {
                    case pca::matrixpcaconst<int32_t>::QVEC : 
                      qvec_rphi = *it;
                      hwmanygot++;
                      hwmanygotrphi++;
                      break;
                    case pca::matrixpcaconst<int32_t>::KVEC :
                      kvec_rphi = *it;
                      hwmanygot++;
                      hwmanygotrphi++;
                      break;
                    case pca::matrixpcaconst<int32_t>::CMTX :
                      cmtx_rphi = *it;
                      hwmanygot++;
                      hwmanygotrphi++;
                      break;
                    case pca::matrixpcaconst<int32_t>::AMTX :
                      amtx_rphi = *it;
                      hwmanygot++;
                      hwmanygotrphi++;
                      break;
                    default:
                      break;
                  }
                }
              } 
            }
          }
        }
      }
    }
  
    if (hwmanygot == 8)
      return true;
    else
    {
      std::cerr << "Found " << hwmanygot << " const instead of 8" << std::endl;
      std::cerr << layersid << " and " << pslayersid << std::endl;
      std::cerr << "charge: " << chargesignin << " eta: " << eta << " pt: " << pt << std::endl;
      std::cerr << "Foundrz " << hwmanygotrz << std::endl;
      std::cerr << "Foundrphi " << hwmanygotrphi << std::endl;
      return false;
    }
  
    // TODO add consistency check for dims
  
    return false;
  }


  bool import_pca_const (
      std::vector<pca::matrixpcaconst<double> > & vct,
      pca::matrixpcaconst<double> & cmtx_rz, 
      pca::matrixpcaconst<double> & qvec_rz, 
      pca::matrixpcaconst<double> & amtx_rz, 
      pca::matrixpcaconst<double> & kvec_rz, 
      pca::matrixpcaconst<double> & cmtx_rphi, 
      pca::matrixpcaconst<double> & qvec_rphi, 
      pca::matrixpcaconst<double> & amtx_rphi, 
      pca::matrixpcaconst<double> & kvec_rphi, 
      double eta, double pt, 
      int chargesignin,
      const std::string & layersid,
      const std::string & pslayersid,
      int towerid)
  {
    int hwmanygot = 0, hwmanygotrphi = 0, hwmanygotrz = 0;;
    std::vector<pca::matrixpcaconst<double> >::const_iterator it = 
      vct.begin();
    for (; it != vct.end(); ++it)
    {
      double ptmin, ptmax, etamin, etamax;
      std::string actuallayids;
      int chargesign;
  
      it->get_ptrange(ptmin, ptmax);
      it->get_etarange(etamin, etamax);
      chargesign = it->get_chargesign();
      actuallayids = it->get_layersids();

      if (towerid == it->get_towerid())
      {
        if (it->get_ttype() == pca::matrixpcaconst<double>::FLOATPT)
        {
          if (it->get_plane_type() == pca::matrixpcaconst<double>::RZ)
          {
            if (actuallayids == pslayersid)
            {
              if ((eta >= etamin) && (eta <= etamax)) 
              {
                switch(it->get_const_type())
                {
                  case pca::matrixpcaconst<double>::QVEC :
                    qvec_rz = *it;
                    hwmanygot++;
                    hwmanygotrz++;
                    break;
                  case pca::matrixpcaconst<double>::KVEC :
                    kvec_rz = *it;
                    hwmanygot++;
                    hwmanygotrz++;
                    break;
                  case pca::matrixpcaconst<double>::CMTX :
                    cmtx_rz = *it;
                    hwmanygot++;
                    hwmanygotrz++;
                    break;
                  case pca::matrixpcaconst<double>::AMTX :
                    amtx_rz = *it;
                    hwmanygot++;
                    hwmanygotrz++;
                    break;
                  default:
                    break;
                }
              } 
            }
          }
          else if (it->get_plane_type() == pca::matrixpcaconst<double>::RPHI)
          {
            if (actuallayids == layersid)
            {
              if (chargesignin == chargesign)
              {
                if ((pt >= ptmin) && (pt <= ptmax))
                {
                  switch(it->get_const_type())
                  {
                    case pca::matrixpcaconst<double>::QVEC : 
                      qvec_rphi = *it;
                      hwmanygot++;
                      hwmanygotrphi++;
                      break;
                    case pca::matrixpcaconst<double>::KVEC :
                      kvec_rphi = *it;
                      hwmanygot++;
                      hwmanygotrphi++;
                      break;
                    case pca::matrixpcaconst<double>::CMTX :
                      cmtx_rphi = *it;
                      hwmanygot++;
                      hwmanygotrphi++;
                      break;
                    case pca::matrixpcaconst<double>::AMTX :
                      amtx_rphi = *it;
                      hwmanygot++;
                      hwmanygotrphi++;
                      break;
                    default:
                      break;
                  }
                }
              } 
            }
          }
        }
      }
    }
  
    if (hwmanygot == 8)
      return true;
    else
    {
      std::cerr << "Found " << hwmanygot << " const instead of 8" << std::endl;
      std::cerr << layersid << " and " << pslayersid << std::endl;
      std::cerr << "charge: " << chargesignin << " eta: " << eta << " pt: " << pt << std::endl;
      std::cerr << "Foundrz " << hwmanygotrz << std::endl;
      std::cerr << "Foundrphi " << hwmanygotrphi << std::endl;
      return false;
    }
  
    // TODO add consistency check for dims
  
    return false;
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
  if (!pca::read_pcacosnt_from_file (pcacontvct_float_, in.c_str()))
  {
    std::cerr << "Error while reading constant from " << in << " read only " << 
      pcacontvct_float_.size() << std::endl;
    return;
  }

  std::vector<pca::matrixpcaconst<double> >::const_iterator it = 
      pcacontvct_float_.begin();
  for (; it != pcacontvct_float_.end(); ++it)
  {
    if (it->get_ttype() != pca::matrixpcaconst<double>::FLOATPT)
    {
      std::cerr << "Wrong PCAconst type " << std::endl;
      return;
    }
  }
  std::cout << "Done" << std::endl;
 
}

void PCATrackFitter::read_integegr_const_filename (const std::string & in)
{
  std::cout << "Reading " << in << std::endl;
  if (!pca::read_pcacosnt_from_file (pcacontvct_integer_, in.c_str()))
  {
    std::cerr << "Error while reading constant from " << in << std::endl;
    return;
  }

  std::vector<pca::matrixpcaconst<int32_t> >::const_iterator it = 
      pcacontvct_integer_.begin();
  for (; it != pcacontvct_integer_.end(); ++it)
  {
    if (it->get_ttype() != pca::matrixpcaconst<int32_t>::INTEGPT)
    {
      std::cerr << "Wrong PCAconst type " << std::endl;
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

  if (hits.size() == 6)
  {
    pca::matrixpcaconst<int32_t> zrv(1, 12), phirv(1, 12);
    std::string layersid, pslayersid;

    if (hits_to_zrpmatrix_integer (ci, si, hits, zrv, phirv, 
          layersid, pslayersid, tow))
    {
      int charge = +1;
      if (track_->getCharge() < 0.0)
       charge = -1;

      // Check the charge TODO
      //charge = -1 * charge;

      double pt_est = track_->getCurve();
      double eta_est = track_->getEta0();
      //double z0_est = track_->getZ0();
      //double phi_est = track_->getPhi0();
      
      pca::matrixpcaconst<int32_t> cmtx_rz(0, 0);
      pca::matrixpcaconst<int32_t> qvec_rz(0, 0); 
      pca::matrixpcaconst<int32_t> amtx_rz(0, 0); 
      pca::matrixpcaconst<int32_t> kvec_rz(0, 0); 
      pca::matrixpcaconst<int32_t> cmtx_rphi(0, 0); 
      pca::matrixpcaconst<int32_t> qvec_rphi(0, 0); 
      pca::matrixpcaconst<int32_t> amtx_rphi(0, 0); 
      pca::matrixpcaconst<int32_t> kvec_rphi(0, 0); 
      
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
                            tow))
      {
        std::cout << "CMTX RZ: " << std::endl;
        dump_element(cmtx_rz, std::cout);
      
        std::cout << "QVEC RZ: " << std::endl;
        dump_element(qvec_rz, std::cout);
      
        std::cout << "CMTX RPHI: " << std::endl;
        dump_element(cmtx_rphi, std::cout);
      
        std::cout << "QVEC RPHI: " << std::endl;
        dump_element(qvec_rphi, std::cout);

        std::cout << "AMTX RZ: " << std::endl;
        dump_element(amtx_rz, std::cout);
      
        std::cout << "KVEC RZ: " << std::endl;
        dump_element(kvec_rz, std::cout);
      
        std::cout << "AMTX RPHI: " << std::endl;
        dump_element(amtx_rphi, std::cout);
      
        std::cout << "KVEC RPHI: " << std::endl;
        dump_element(kvec_rphi, std::cout);
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

  //std::cout << "PCA::fit tow: " << tow << " hits size: " << 
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
                            tow))
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
                            tow))
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

