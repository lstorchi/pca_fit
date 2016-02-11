/*
C++ implementation of the PCA fitter

L.Storchi, A.Modak, S.R.Chowdhury : 2016
*/
#include "../interface/PCATrackFitter.h"

// to store PCA constants very quick, quite dirty and very ugly (may be a class or ....)
// is there a 2d array class in cmssw ? 
#define PCA_RZ_RDIM 2
#define PCA_RZ_CDIM 6
#define PCA_RPHI_RDIM 2
#define PCA_RPHI_CDIM 12
namespace 
{
  double c_18_rz_1[PCA_RZ_RDIM][PCA_RZ_CDIM] = {};
  double q_18_rz_1[PCA_RZ_RDIM] = {};
  double c_18_rphi_1_cgt0[PCA_RPHI_RDIM][PCA_RPHI_CDIM] = {};
  double q_18_rphi_1_cgt0[PCA_RPHI_RDIM] = {};

  void get_c_matrix_rz (double eta, int towerid, 
      double c[PCA_RZ_RDIM][PCA_RZ_CDIM])
  {
    if (towerid == 18)
    {
      if ((eta >= -0.6) && (eta <= -0.55))
      {
        // TODO read and store
        for (int i=0; i<PCA_RZ_RDIM; ++i)
          for (int j=0; j<PCA_RZ_CDIM; ++j)
            c[i][j] = c_18_rz_1[i][j];
      }
      else 
      {
        // TODO
      }
    }
    else
    {
      // TODO
    }
  }

  void get_c_matrix_rphi (double pt, int charge, 
      int towerid, double c[PCA_RPHI_RDIM][PCA_RPHI_CDIM])
  {
    if (charge > 0)
    {
      if (towerid == 18)
      {
        if ((pt >= 3.0) && (pt <= 7))
        {
          // TODO read and store
          for (int i=0; i<PCA_RPHI_RDIM; ++i)
            for (int j=0; j<PCA_RPHI_CDIM; ++j)
              c[i][j] = c_18_rphi_1_cgt0[i][j];
        }
        else 
        {
          // TODO
        }
      }
      else
      {
        // TODO
      }
    }
    else
    {
      // TODO
    }
  }


  void get_q_vector_rz (double eta, int towerid, double q[PCA_RZ_RDIM])
  {
    if (towerid == 18)
    {
      if ((eta >= -0.6) && (eta <= -0.55))
      {
        // TODO read and store
        for (int i=0; i<PCA_RZ_RDIM; ++i)
          q[i] = q_18_rz_1[i];
      }
      else 
      {
        // TODO
      }
    }
    else
    {
      // TODO
    }
  }

  void get_q_vector_rphi (double pt, int charge,  
      int towerid, double q[PCA_RPHI_RDIM])
  {
    if (charge > 0)
    {
      if (towerid == 18)
      {
        if ((pt >= 3.0) && (pt <= 7))
        {
          // TODO read and store
          for (int i=0; i<PCA_RPHI_RDIM; ++i)
            q[i] = q_18_rphi_1_cgt0[i];
        }
        else 
        {
          // TODO
        }
      }
      else
      {
        // TODO
      }
    }
    else 
    {
      // TODO
    }
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

  std::cout << "In PCA::fit tow: " << tow << std::endl;

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

  std::cout << "tracks_.size() : " << tracks_.size() << std::endl;

  for(unsigned int tt=0; tt<tracks_.size(); ++tt)
  {
    //vector<int> stubids = tracks_[tt]->getStubs(); // TODO are they the hits idx ? 

    std::cout << "hits.size() : " << hits.size() << std::endl;

    if (hits.size() == 6)
    {
      std::vector<double> zrv, phirv;
      int charge = tracks_[tt]->getCharge(); 

      for(unsigned int idx = 0; idx < hits.size(); ++idx)
      {
        double xi =  hits[idx]->getX()*ci+ hits[idx]->getY()*si;
        double yi = -hits[idx]->getX()*si+ hits[idx]->getY()*ci;

        double zi = hits[idx]->getZ();
        double ri = sqrt(xi*xi+yi*yi);
        double pi = atan2(yi,xi);

        std::cout << (int) hits[idx]->getLayer() << std::endl;
        std::cout << xi << " " << yi << " " << zi << std::endl;
        std::cout << ri << " " << pi << std::endl;

        zrv.push_back(zi);
        zrv.push_back(ri);

        phirv.push_back(pi);
        phirv.push_back(ri);
      }

      double pt_est, eta_est;

      pt_est = tracks_[tt]->getCurve();
      eta_est = tracks_[tt]->getEta0();

      double c_rz[PCA_RZ_RDIM][PCA_RZ_CDIM];
      double q_rz[PCA_RZ_RDIM];

      double c_rphi[PCA_RPHI_RDIM][PCA_RPHI_CDIM];
      double q_rphi[PCA_RPHI_RDIM];

      double cottheta = 0.0; // eta
      double z0 = 0.0;
      get_c_matrix_rz (eta_est, tow, c_rz);
      get_q_vector_rz (eta_est, tow, q_rz);

      // TODO checkit
      cottheta = q_rz[0];
      z0 = q_rz[1];
      for (int i=0; i<PCA_RPHI_CDIM; ++i)
      {
        cottheta += c_rz[0][i]*zrv[i];
        z0 += c_rz[1][i]*zrv[i];
      }
 
      double coverpt = 0.0; // pt
      double phi = 0.0;
      get_c_matrix_rphi (pt_est, charge, tow, c_rphi);
      get_q_vector_rphi (pt_est, charge, tow, q_rphi);

      // TODO checkit
      coverpt = q_rphi[0];
      phi = q_rphi[1];
      for (int i=0; i<PCA_RPHI_CDIM; ++i)
      {
        coverpt += c_rphi[0][i]*phirv[i];
        z0 += c_rphi[1][i]*phirv[i];
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
      // TODO non 6 layers 
    }
  }
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
