#include "PrincipalFitGenerator.h"

#include "PatternFinder.h"

PrincipalFitGenerator::PrincipalFitGenerator(string f, SectorTree *s){
  inputDirectory = f;
  st = s;
}

TChain* PrincipalFitGenerator::createTChain(){

  cout<<"Loading files from "<<inputDirectory<<" (this may take several minutes)..."<<endl;

  TChain* TT = new TChain("BankStubs");

  TSystemDirectory dir(inputDirectory.c_str(), inputDirectory.c_str());
  TList *files = dir.GetListOfFiles();

  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);

    while((file=(TSystemFile*)next())){
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(".root")) {
	TT->Add((inputDirectory + fname.Data()).c_str());
      }
    }
  }

  p_m_stub_modid = &m_stub_modid; 
  p_m_stub_ptGEN = &m_stub_ptGEN;  
  p_m_stub_phiGEN = &m_stub_phiGEN;  
  p_m_stub_etaGEN = &m_stub_etaGEN;  
  p_m_stub_x = &m_stub_x;
  p_m_stub_y = &m_stub_y;
  p_m_stub_z = &m_stub_z;
  p_m_stub_pdg = &m_stub_pdg;
  p_m_stub_x0 = &m_stub_x0;
  p_m_stub_y0 = &m_stub_y0;
  p_m_stub_z0 = &m_stub_z0;

  TT->SetBranchAddress("STUB_n",         &m_stub);
  TT->SetBranchAddress("STUB_modid",     &p_m_stub_modid);
  TT->SetBranchAddress("STUB_ptGEN",     &p_m_stub_ptGEN);
  TT->SetBranchAddress("STUB_phiGEN",    &p_m_stub_phiGEN);
  TT->SetBranchAddress("STUB_etaGEN",    &p_m_stub_etaGEN);
  TT->SetBranchAddress("STUB_x",         &p_m_stub_x);
  TT->SetBranchAddress("STUB_y",         &p_m_stub_y);
  TT->SetBranchAddress("STUB_z",         &p_m_stub_z);
  TT->SetBranchAddress("STUB_pdgGEN",     &p_m_stub_pdg);
  TT->SetBranchAddress("STUB_xGEN",    &p_m_stub_x0);
  TT->SetBranchAddress("STUB_yGEN",    &p_m_stub_y0);
  TT->SetBranchAddress("STUB_zGEN",    &p_m_stub_z0);

  int nb_entries = TT->GetEntries();
  cout<<nb_entries<<" events found."<<endl;

  return TT;
}


// covariance of coodinates ?

void PrincipalFitGenerator::generatePrincipal(map<int,pair<float,float> > eta_limits, float min_pt, float max_pt, float min_eta, float max_eta){

  if(st->getNbSectors()==0){
    cout<<"No sector found"<<endl;
    return;
  }

  map< int, vector<int> > detector_config = Sector::readConfig("detector.cfg");
  //Get the layers IDs
  vector<int> tracker_layers;
  Sector* first_sector = st->getAllSectors()[0];
  for(int i=0;i<first_sector->getNbLayers();i++){
    tracker_layers.push_back(first_sector->getLayerID(i));
  }
  

int nb_ladders = -1;
  if(detector_config.find(first_sector->getLayerID(0))!=detector_config.end()){
    nb_ladders = detector_config[first_sector->getLayerID(0)][1];
  }
  else{
    cout<<"We do not know the number of ladders in layer "<<first_sector->getLayerID(0)<<endl;
    return;
  }

  TChain* TT = createTChain();
  int n_entries_TT = TT->GetEntries();

  int evtIndex=0;

  int layers[tracker_layers.size()];
  vector<int> ladder_per_layer(tracker_layers.size());
  vector<int> module_per_layer(tracker_layers.size());

  int max_tracks_number = 10000000;

  while(evtIndex<n_entries_TT){
    TT->GetEntry(evtIndex++);

    if(evtIndex%100000==0)
      cout<<"PCA : Event "<<evtIndex<<endl;

    if(!selectTrack(tracker_layers, layers, ladder_per_layer, module_per_layer, eta_limits, min_pt, max_pt, min_eta, max_eta))
      continue;

    Sector* sector = st->getSector(ladder_per_layer, module_per_layer);
    if(sector==NULL){
      //cout<<"No sector found"<<endl;
      continue;
    }

    int modid = m_stub_modid[layers[0]];
    int first_layer = modid/1000000;
    modid = modid-first_layer*1000000;
    int first_ladder = modid/10000;

    double sec_phi = (2*M_PI/nb_ladders)*(sector->getLadderCode(tracker_layers[0],CMSPatternLayer::getLadderCode(tracker_layers[0],first_ladder)));

    PrincipalTrackFitter* fitter = (PrincipalTrackFitter*)sector->getFitter();
    if(fitter==NULL){
      cout<<"No fitter associated to the sector, creating a default one!"<<endl;
      fitter = new PrincipalTrackFitter(sector->getNbLayers(), 1000);
      sector->setFitter(fitter);
    }
    fitter->setPhiRotation(sec_phi);

    double data_coord[3*tracker_layers.size()];
    int data_tracker[3*tracker_layers.size()];

    for(unsigned int i=0;i<tracker_layers.size();i++){
      int j = layers[i];
      
      int value = m_stub_modid[j];
      int layer = value/1000000;
      value = value-layer*1000000;
      int ladder = value/10000;
      value = value-ladder*10000;
      int module = value/100;

      data_coord[i*3]=m_stub_x[j];
      data_coord[i*3+1]=m_stub_y[j];
      data_coord[i*3+2]=m_stub_z[j];

      data_tracker[i*3]=sector->getLayerIndex(layer);
      data_tracker[i*3+1]=sector->getLadderCode(layer,ladder);
      data_tracker[i*3+2] = sector->getModuleCode(layer, CMSPatternLayer::getLadderCode(layer,ladder), CMSPatternLayer::getModuleCode(layer,module));
    }

    fitter->addTrackForPrincipal(data_tracker, data_coord);
    if(fitter->hasPrincipalParams())//The params are computed->we stop the process
      break;

    if(evtIndex>max_tracks_number){
      fitter->forcePrincipalParamsComputing();
    }

  }

  delete TT;
}

// covariance of coodinates vs parameters ?

void PrincipalFitGenerator::generateMultiDim(map<int,pair<float,float> > eta_limits, float min_pt, float max_pt, float min_eta, float max_eta){

  if(st->getNbSectors()==0){
    cout<<"No sector found"<<endl;
    return;
  }

  map< int, vector<int> > detector_config = Sector::readConfig("detector.cfg");

  //Get the layers IDs
  vector<int> tracker_layers;
  Sector* first_sector = st->getAllSectors()[0];
  cout<<"We use layers : ";
  for(int i=0;i<first_sector->getNbLayers();i++){
    tracker_layers.push_back(first_sector->getLayerID(i));
    cout<<first_sector->getLayerID(i)<<",";
  }
  cout<<endl;

  /*
  int nb_ladders = -1;
  if(detector_config.find(first_sector->getLayerID(0))!=detector_config.end()){
    nb_ladders = detector_config[first_sector->getLayerID(0)][1];
  }
  else{
    cout<<"We do not know the number of ladders in layer "<<first_sector->getLayerID(0)<<endl;
    return;
  }
  */

  TChain* TT = createTChain();
  int n_entries_TT = TT->GetEntries();

  int evtIndex=0;

  int layers[tracker_layers.size()];
  vector<int> ladder_per_layer(tracker_layers.size());
  vector<int> module_per_layer(tracker_layers.size());

  int max_tracks_number = 10000000;

  cerr << "n_entries_TT: " << n_entries_TT << std::endl;
  
  while(evtIndex<n_entries_TT){
    TT->GetEntry(evtIndex++);

    if(evtIndex%100000==0)
      cout<<"MultiDim : Event "<<evtIndex<<endl;

    if(!selectTrack(tracker_layers, layers, ladder_per_layer, module_per_layer, eta_limits, min_pt, max_pt, min_eta, max_eta))
      continue;
    
    //cout<<"track is complete, searching for a sector"<<endl;

    Sector* sector = st->getSector(ladder_per_layer, module_per_layer);
    if(sector==NULL){
      //cout<<"No sector found"<<endl;
      continue;
    }

    PrincipalTrackFitter* fitter = (PrincipalTrackFitter*)sector->getFitter();
    if(fitter==NULL){
      cout<<"No fitter associated to the sector!"<<endl;
      break;
    }
    if(!fitter->hasPrincipalParams()){
      cout<<"No principal parameters found!"<<endl;
      break;
    }

    double data_coord[3*tracker_layers.size()];
    int data_tracker[3*tracker_layers.size()];
  
    cerr << evtIndex << " " << tracker_layers.size() <<  std::endl;

    for(unsigned int i=0;i<tracker_layers.size();i++){
      int j = layers[i];
      
      int value = m_stub_modid[j];
      int layer = value/1000000;
      value = value-layer*1000000;
      int ladder = value/10000;
      value = value-ladder*10000;
      int module = value/100;

      //data_coord[i*3]=m_stub_x[j]*cos(sec_phi)+m_stub_y[j]*sin(sec_phi);
      //data_coord[i*3+1]=-m_stub_x[j]*sin(sec_phi)+m_stub_y[j]*cos(sec_phi);
      data_coord[i*3]=m_stub_x[j];
      data_coord[i*3+1]=m_stub_y[j];
      data_coord[i*3+2]=m_stub_z[j];

      data_tracker[i*3]=sector->getLayerIndex(layer);
      data_tracker[i*3+1]=sector->getLadderCode(layer,CMSPatternLayer::getLadderCode(layer,ladder));
      data_tracker[i*3+2] = sector->getModuleCode(layer, CMSPatternLayer::getLadderCode(layer,ladder), CMSPatternLayer::getModuleCode(layer,module));

      cerr << data_coord[i*3] << " " << data_coord[i*3+1] << " " << data_coord[i*3+2] << " " <<
	 data_tracker[i*3] << " " << data_tracker[i*3+1] << " " << data_tracker[i*3+2] << endl;
    }

    int charge = m_stub_pdg[layers[0]]/abs(m_stub_pdg[layers[0]]);
    double d0 = (m_stub_y0[layers[0]]-tan(m_stub_phiGEN[layers[0]])*m_stub_x0[layers[0]])*cos(m_stub_phiGEN[layers[0]]);

    double values[5];
    values[0] = charge*m_stub_ptGEN[layers[0]];//PT
    values[1] = m_stub_phiGEN[layers[0]];
    values[2] = d0;
    values[3] = m_stub_etaGEN[layers[0]];
    values[4] = m_stub_z0[layers[0]];

    cerr << values[0] << " " << values[1] << " " << values[2] << " " 
         << values[3] << " " << values[4] << endl;

    fitter->addTrackForMultiDimFit(data_tracker, data_coord, values);

    if(fitter->hasMultiDimFitParams())//The params are computed->we stop the process
      break;

    if(evtIndex>max_tracks_number){
      fitter->forceMultiDimFitParamsComputing();
    }
  }

  delete TT;
}

bool PrincipalFitGenerator::selectTrack(vector<int> &tracker_layers, int* layers, vector<int> &ladder_per_layer, vector<int> &module_per_layer,
					map<int,pair<float,float> > &eta_limits, float min_pt, float max_pt, float min_eta, float max_eta){

   //initialize arrays
    for(unsigned int j=0;j<tracker_layers.size();j++){
      layers[j]=-1;
    }
    for(unsigned int j=0;j<tracker_layers.size();j++){
      ladder_per_layer[j]=-1;
    }
    for(unsigned int j=0;j<tracker_layers.size();j++){
      module_per_layer[j]=-1;
    }

    float current_eta = -1;

    for(int j=0;j<m_stub;j++){

      if(m_stub_etaGEN[j]<min_eta){// eta of the generating particule is bellow the threshold -> we do not use it
      	return false;
      }
      if(m_stub_etaGEN[j]>max_eta){// eta of the generating particule is above the threshold -> we do not use it
      	return false;
      }

      if(m_stub_ptGEN[j]<min_pt){// The PT of the generating particule is below the minimum required -> we do not use it
	return false;
      }
      if(m_stub_ptGEN[j]>max_pt){// The PT of the generating particule is above the maximum accepted -> we do not use it
	return false;
      }

      int value = m_stub_modid[j];
      int modid_layer = value/1000000;
      value = value-modid_layer*1000000;
      int modid_ladder = value/10000;
      value = value-modid_ladder*10000;
      int modid_module = value/100;

      int layer_position=-1;
      for(unsigned int cpt=0;cpt<tracker_layers.size();cpt++){
	if(modid_layer==tracker_layers[cpt]){
	  layer_position=(int)cpt;
	  break;
	}
      }
      
      if(layer_position!=-1){ // is this layer in the layer list?
	layers[layer_position]=j;
	ladder_per_layer[layer_position]=CMSPatternLayer::getLadderCode(modid_layer,modid_ladder);
	short module = -1;
	module = CMSPatternLayer::getModuleCode(modid_layer, modid_module);
	module_per_layer[layer_position]=module;
      }
      
      current_eta = m_stub_etaGEN[j];
    }

    /**************************************
    Selection on the stubs/layer
    We need at least one stub per layer
    **************************************/
    bool missing_stub = false;
    for(unsigned int j=0;j<tracker_layers.size();j++){
      if(layers[j]==-1){
	missing_stub=true;
	if(st->getNbSectors()==1){
	  if(eta_limits.find(tracker_layers[j])!=eta_limits.end()){//we have eta boundaries for this layer
	    pair<float,float> limits = eta_limits[tracker_layers[j]];
	    if(current_eta<limits.first || current_eta>limits.second){ // we are outside the eta limits for this layer
	      //cout<<"missing hit on layer "<<tracker_layers[j]<<" for track with eta="<<current_eta<<endl;
	      layers[j]=-2;
	      //we put a ladder and a module just to be inside the sector
	      ladder_per_layer[j]=st->getAllSectors()[0]->getLadders(j)[0];
	      module_per_layer[j]=st->getAllSectors()[0]->getModules(j,ladder_per_layer[j])[0];
	      //cout<<"Add stub for sector : "<<ladder_per_layer[j]<<" / "<<module_per_layer[j]<<endl;
	      //debug=true;
	      missing_stub=false;//we will create a fake superstrip, so the stub is not missing
	    }
	  }
	}
      }
      if(missing_stub)
	break;
    }

    if(missing_stub){
      return false;//no stub on each layer -> drop the event    
    }

    return true;
}

void PrincipalFitGenerator::generate(map<int,pair<float,float> > eta_limits, float min_pt, float max_pt, float min_eta, float max_eta){
  generatePrincipal(eta_limits, min_pt, max_pt, min_eta, max_eta);
  generateMultiDim(eta_limits, min_pt, max_pt, min_eta, max_eta);
}
