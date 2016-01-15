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

void PCATrackFitter::setTracks(std::vector<Track*>)
{

}

void PCATrackFitter::fit(vector<Hit*> hits)
{
  // TODO
}

void PCATrackFitter::fit(){

  vector<Hit*> activatedHits;

  //////// Get the list of unique stubs from the tracks ///////////
  // TODO

  fit(activatedHits);
}

TrackFitter* PCATrackFitter::clone()
{
  PCATrackFitter* fit = new PCATrackFitter(nb_layers);
  return fit;
}
