#include "../interface/PCATrackFitter.h"

//
// Use the unnamed namespace for things only used in this file
//
namespace
{
  /** @brief A guard object that will delete all of the objects in a container when
   * it goes out of scope.
   *
   * This should have the same lifetime as the vector who's elements it will delete.
   * The template parameter should be the container type, not the contained object.
   *
   * This is so that a vector of 'new'ed pointers can have each element deleted when
   * a function returns, without having to worry about where it returns.
   * It would have been much easier to use an edm::OwnVector but the vector is passed
   * around and it would mean changing the class interface. This seemed to be the least
   * invasive way to achieve the required result.
   */
  template<class T>
  class DeleteContainerElementsGuard
  {
    public:
      DeleteContainerElementsGuard( T& container ) : 
        container_(container) 
      {
      }

      ~DeleteContainerElementsGuard() 
      { 
        for( auto pElement : container_ ) 
          delete pElement; 
      }

    protected:
      T& container_;
  };
}

PCATrackFitter::PCATrackFitter():TrackFitter(0)
{
  verboseLevel_ = 0;
  event_counter_ = 0;
  road_id_ = 0;

  initialize();
}

PCATrackFitter::PCATrackFitter(int nb):TrackFitter(nb)
{
  verboseLevel_ = 0;
  event_counter_ = 0;
  road_id_ = 0;

  nb_layers = nb;

  initialize();
}

PCATrackFitter::~PCATrackFitter()
{
}

void PCATrackFitter::fit(vector<Hit*> hits_)
{
  if ( hits_.size() != (unsigned int) nb_layers )
  {
    cout << "*** ERROR in PCATrackFitter::fit() at event/road = " << event_counter_ << "/" 
	 << road_id_ << ": only " << hits_.size() << " stubs found, fit aborted!" << endl;
    return;
  }

  // TODO

}

void PCATrackFitter::fit()
{

  vector<Hit*> activatedHits;

  //////// Get the list of unique stubs from all the patterns ///////////
  set<int> ids;
  int total=0;
  
  for(unsigned int i=0;i<patterns.size();i++)
  {
    vector<Hit*> allHits = patterns[i]->getHits();
    total+=allHits.size();
    
    for(unsigned int j=0;j<allHits.size();j++)
    {
      pair<set<int>::iterator,bool> result = ids.insert(allHits[j]->getID());
      if(result.second==true)
	activatedHits.push_back(allHits[j]);
    }
  }

  fit(activatedHits);
}


void PCATrackFitter::mergePatterns()
{
  //cout<<"Merging of patterns not implemented"<<endl;
}


void PCATrackFitter::mergeTracks()
{
}

TrackFitter* PCATrackFitter::clone()
{
  PCATrackFitter* fit = new PCATrackFitter(nb_layers);
  fit->setPhiRotation(sec_phi);
  fit->setSectorID(sector_id);

  return fit;
}

void PCATrackFitter::initialize()
{

}
