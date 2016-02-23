/*
C++ implementation of the PCA fitter

L.Storchi: 2016
*/
#include "../interface/PCATrackFitter.h"
#include "../interface/pcaconst.hpp"

namespace 
{
  #define LAYIDDIM 6
  int stdlayersid[LAYIDDIM] = {5, 6, 7, 8, 9, 10};

  /* quick and very dirty */
  void dump_element (const pca::matrixpcaconst<double> & in, std::ostream & out)
  {
    out.precision(6);
    for (unsigned int i = 0; i<in.n_rows(); ++i)
    {
      for (unsigned int j = 0; j<in.n_cols(); ++j)
        out << std::scientific << in(i, j) << " ";
      out << std::endl;
    }
  }

  bool hits_to_zrpmatrix (double ci, double si, 
      const std::vector<Hit*> & hits, 
      const std::vector<int> & stubs,
      pca::matrixpcaconst<double> & zrv, 
      pca::matrixpcaconst<double> & phirv)
  {
    std::vector<int> layerids;
    std::vector<double> riv, piv, ziv;

    /*
    for(unsigned int idx = 0; idx < hits.size(); ++idx)
    {
    */
    for (unsigned int sti=0; sti<stubs.size(); sti++) 
    {
      unsigned int idx = stubs[sti] - 1;
      // TODO double check this
      double xi =  hits[idx]->getX()*ci+ hits[idx]->getY()*si;
      double yi = -hits[idx]->getX()*si+ hits[idx]->getY()*ci;

      xi = hits[idx]->getX();
      yi = hits[idx]->getY();

      double zi = hits[idx]->getZ();
      double ri = sqrt(xi*xi+yi*yi);
      double pi = atan2(yi,xi);

      /*
      std::cout << "Layer  " << (int) hits[idx]->getLayer() << std::endl;
      std::cout << "PCAxyz " << xi << " " << yi << " " << zi << std::endl;
      std::cout << "PCArpz " << ri << " " << pi << " " << zi << std::endl;
      */

      layerids.push_back((int)hits[idx]->getLayer());

      riv.push_back(ri);
      piv.push_back(pi);
      ziv.push_back(zi);
    }

    double z1 = 0.0, z3 = 0.0, r1 = 0.0, r3 = 0.0;

    int counter = 0;
    for (int i=0; i<LAYIDDIM; ++i)
    {
      // readoder hits quick test
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
    double intercept = z1 - slope*r1;

    std::cout << "Intercept est z0: " << intercept << std::endl;
 
    return true;
  }

  bool import_pca_const (const std::string & cfname, 
      pca::matrixpcaconst<double> & cmtx_rz, 
      pca::matrixpcaconst<double> & qvec_rz, 
      pca::matrixpcaconst<double> & amtx_rz, 
      pca::matrixpcaconst<double> & kvec_rz, 
      pca::matrixpcaconst<double> & cmtx_rphi, 
      pca::matrixpcaconst<double> & qvec_rphi, 
      pca::matrixpcaconst<double> & amtx_rphi, 
      pca::matrixpcaconst<double> & kvec_rphi, 
      double eta, double pt, 
      int chargesignin)
  {
    std::vector<pca::matrixpcaconst<double> > vct;
    if (pca::read_pcacosnt_from_file (vct, cfname.c_str()))
    {
      int hwmanygot = 0;
      std::vector<pca::matrixpcaconst<double> >::const_iterator it = 
        vct.begin();
      for (; it != vct.end(); ++it)
      {
        double ptmin, ptmax, etamin, etamax;
        int chargesign;
  
        it->get_ptrange(ptmin, ptmax);
        it->get_etarange(etamin, etamax);
        chargesign = it->get_chargesign();
  
        if (it->get_plane_type() == pca::matrixpcaconst<double>::RZ)
        {
          if ((eta >= etamin) && (eta <= etamax)) 
          {
            switch(it->get_const_type())
            {
              case pca::matrixpcaconst<double>::QVEC :
                qvec_rz = *it;
                hwmanygot++;
                break;
              case pca::matrixpcaconst<double>::KVEC :
                kvec_rz = *it;
                hwmanygot++;
                break;
              case pca::matrixpcaconst<double>::CMTX :
                cmtx_rz = *it;
                hwmanygot++;
                break;
              case pca::matrixpcaconst<double>::AMTX :
                amtx_rz = *it;
                hwmanygot++;
                break;
              default:
                break;
            }
          } 
        }
        else if (it->get_plane_type() == pca::matrixpcaconst<double>::RPHI)
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
                  break;
                case pca::matrixpcaconst<double>::KVEC :
                  kvec_rphi = *it;
                  hwmanygot++;
                  break;
                case pca::matrixpcaconst<double>::CMTX :
                  cmtx_rphi = *it;
                  hwmanygot++;
                  break;
                case pca::matrixpcaconst<double>::AMTX :
                  amtx_rphi = *it;
                  hwmanygot++;
                  break;
                default:
                  break;
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
        return false;
      }
    }
  
    // TODO add consistency check for dims
  
    return false;
  }
}

PCATrackFitter::PCATrackFitter():TrackFitter(0)
{
}

PCATrackFitter::PCATrackFitter(int nb):TrackFitter(nb)
{
}

PCATrackFitter::~PCATrackFitter()
{
}

void PCATrackFitter::initialize()
{
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

void PCATrackFitter::setTracks(const std::vector<Track*> & intc)
{
  tracks_ = intc;
}

void PCATrackFitter::fit(vector<Hit*> hits)
{
  int tow = sector_id; // The tower ID, necessary to get the phi shift

  std::cout << "In PCA::fit tow: " << tow << " hits size: " << 
    hits.size() << std::endl;
  
  double sec_phi = 0;
  switch (tow%8)
  {
    case 0:
      sec_phi = 0.4;
      break;
    case 1:
      sec_phi = 1.2;
      break;
    case 3:
      sec_phi = 2.0;
      break;
    case 4:
      sec_phi = 2.7;
      break;
    case 5:
      sec_phi = -2.0;
      break;
    case 6:
      sec_phi = -1.2;
      break;
    case 7:
      sec_phi = -0.4;
      break;
  }
  
  double ci = cos(sec_phi);
  double si = sin(sec_phi);

  //std::cout << "tracks_.size() : " << tracks_.size() << std::endl;
  
  /*
  std::ofstream xyzfile;
  xyzfile.open ("PCAxyzfile.txt");

  std::ofstream rpzfile;
  rpzfile.open ("PCArphizfile.txt");
  */

  std::cout << "Track size: " << tracks_.size() << std::endl;

  for(unsigned int tt=0; tt<tracks_.size(); ++tt)
  {
    //std::cout << "hits.size() : " << hits.size() << std::endl;
    std::vector<int> stubs = tracks_[tt]->getStubs(); // TODO are they the hits idx ? 

    if (stubs.size() == 6)
    {
      pca::matrixpcaconst<double> zrv(1, 12), phirv(1, 12);

      if (hits_to_zrpmatrix (ci, si, hits, stubs, zrv, phirv))
      {
        int charge;
        double pt_est, eta_est, z0_est, phi_est;
        charge = tracks_[tt]->getCharge(); 
        pt_est = tracks_[tt]->getCurve();
        eta_est = tracks_[tt]->getEta0();
        z0_est = tracks_[tt]->getZ0();
        phi_est = tracks_[tt]->getPhi0();
        
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
        
        std::string cfname = "./barrel_tow18_pca_const.txt";
        
        pca::matrixpcaconst<double> cmtx_rz(0, 0);
        pca::matrixpcaconst<double> qvec_rz(0, 0); 
        pca::matrixpcaconst<double> amtx_rz(0, 0); 
        pca::matrixpcaconst<double> kvec_rz(0, 0); 
        pca::matrixpcaconst<double> cmtx_rphi(0, 0); 
        pca::matrixpcaconst<double> qvec_rphi(0, 0); 
        pca::matrixpcaconst<double> amtx_rphi(0, 0); 
        pca::matrixpcaconst<double> kvec_rphi(0, 0); 
        
        if (import_pca_const (cfname, 
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
                              charge))
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
          
          double pt = (double)(charge)/coverpt;
          
          // TODO: checkit theta to eta 
          double eta = 0.0e0;
          double theta = atan(1.0e0 / cottheta); 
          double tantheta2 = tan (theta/2.0e0); 
          if (tantheta2 < 0.0)
            eta = 1.0e0 * log (-1.0e0 * tantheta2);
          else
            eta = -1.0e0 * log (tantheta2);
          
          Track* fit_track = new Track();
        
          std::cout << " pt : " << pt << " ==> " << pt_est << std::endl;
          std::cout << " phi: " << phi << " ==> " << phi_est << std::endl; 
          std::cout << " eta: " << eta << " ==> " << eta_est << std::endl;
          std::cout << " z0 : " << z0 << " ==> " << z0_est << std::endl;
          
          fit_track->setCurve(pt);
          fit_track->setPhi0(phi);
          fit_track->setEta0(eta);
          fit_track->setZ0(z0);
                          
          for(unsigned int idx = 0; idx < hits.size(); ++idx)
            fit_track->addStubIndex(idx);
          
          tracks.push_back(fit_track);
        }
        else 
        {
          std::cerr << "error while reading PCA const" << std::endl;
          std::cout << "error while reading PCA const" << std::endl;
        }
      } 
      else
      {
        std::cerr << "error in coord conv" << std::endl;
        std::cout << "error in coord conv" << std::endl;
      }

    }
    else
    {
      // TODO non 6 layers 
    }
  }

  /*
  xyzfile.close();
  rpzfile.close();

  std::cout << "Close files" << std::endl;
  */
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

