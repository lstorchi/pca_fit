/*
C++ implementation of the PCA fitter

L.Storchi, A.Modak, S.R.Chowdhury : 2016
*/

#include "../interface/PCATrackFitter.h"

PCATrackFitter::PCATrackFitter():TrackFitter(0){

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

void PCATrackFitter::setTracks(std::vector<Track*> & intc)
{
  tracks_ = intc;
}

void PCATrackFitter::fit(vector<Hit*> hits)
{
  int tow = sector_id; // The tower ID, necessary to get the phi shift

  float sec_phi = 0;
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
  
  float ci = cos(sec_phi);
  float si = sin(sec_phi);

  for(unsigned int tt=0; tt<tracks_.size(); ++tt)
  {
    vector<int> stubids = tracks_[tt]->getStubs();
    if (stubids.size() == 6)
    {
      std::vector<float> zv, rv, pv;

      for(unsigned int idx_i=0; idx_i<stubids.size(); ++idx_i)
      {
        float xi = hits[idx_i]->getX()*ci+ hits[idx_i]->getY()*si;
        float yi = -hits[idx_i]->getX()*si+ hits[idx_i]->getY()*ci;

        float zi = hits[idx_i]->getZ();
        float ri = sqrt(xi*xi+yi*yi);
        float pi = atan2(yi,xi);

        zv.push_back(zi);
        rv.push_back(ri);
        pv.push_back(pi);

        std::cout << xi << yi << zi << ri << pi << std::endl;
      }

      // TODO
    }
    else
    {
      // TODO non 6 layers 
    }
  }
}

void PCATrackFitter::fit(){

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
