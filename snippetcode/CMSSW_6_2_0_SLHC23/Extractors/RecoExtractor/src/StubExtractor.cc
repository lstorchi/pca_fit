#include "../interface/StubExtractor.h"


StubExtractor::StubExtractor(std::string CLUS_tag,std::string CLUS_name,std::string STUB_tag,std::string STUB_name, bool doTree)
{

  m_CLUS_tag   = CLUS_tag;
  m_STUB_tag   = STUB_tag;

  m_CLUS_name  = CLUS_name;
  m_STUB_name  = STUB_name;

  m_OK = false;
  n_tot_evt=0;

  // Tree definition
 
  m_clus_x       = new  std::vector<float>;
  m_clus_y       = new  std::vector<float>;
  m_clus_z       = new  std::vector<float>; 
  m_clus_e       = new  std::vector<float>;
  m_clus_layer   = new  std::vector<int>;
  m_clus_module  = new  std::vector<int>;  
  m_clus_ladder  = new  std::vector<int>;  
  m_clus_seg     = new  std::vector<int>; 
  m_clus_strip   = new  std::vector<float>; 
  m_clus_used    = new  std::vector<int>;  
  m_clus_sat     = new  std::vector<int>; 
  m_clus_nstrips = new  std::vector<int>; 
  m_clus_matched = new  std::vector<int>;  
  m_clus_PS      = new  std::vector<int>;  
  m_clus_nrows   = new  std::vector<int>;
  m_clus_pid     = new  std::vector<int>;  
  m_clus_pdgID   = new  std::vector<int>;  
  m_clus_ptGEN   = new  std::vector<float>; 
  m_clus_rank_PIX= new  std::vector<int>;  
  m_clus_rank_STR= new  std::vector<int>; 
  m_clus_rank_PIX_4= new  std::vector<int>; 
  m_clus_rank_STR_4= new  std::vector<int>; 
  m_clus_rank_PIX_8= new  std::vector<int>; 
  m_clus_rank_STR_8= new  std::vector<int>;  
  m_clus_tp      = new  std::vector<int>;  
  m_clus_pix     = new  std::vector<std::vector<int> >; 
  m_clus_mult    = new  std::vector<std::vector<int> >; 

  m_stub_pt      = new  std::vector<float>;
  m_stub_ptMC    = new  std::vector<float>; 
  m_stub_pxGEN   = new  std::vector<float>; 
  m_stub_pyGEN   = new  std::vector<float>; 
  m_stub_etaGEN  = new  std::vector<float>; 
  m_stub_X0      = new  std::vector<float>; 
  m_stub_Y0      = new  std::vector<float>; 
  m_stub_Z0      = new  std::vector<float>; 
  m_stub_PHI0    = new  std::vector<float>; 
  m_stub_layer   = new  std::vector<int>; 
  m_stub_module  = new  std::vector<int>;  
  m_stub_ladder  = new  std::vector<int>; 
  m_stub_seg     = new  std::vector<int>;  
  m_stub_chip    = new  std::vector<int>;  
  m_stub_strip   = new  std::vector<float>; 
  m_stub_x       = new  std::vector<float>;  
  m_stub_y       = new  std::vector<float>;  
  m_stub_z       = new  std::vector<float>;  
  m_stub_clust1  = new  std::vector<int>;  
  m_stub_clust2  = new  std::vector<int>;  
  m_stub_cw1     = new  std::vector<int>;  
  m_stub_cw2     = new  std::vector<int>;  
  m_stub_deltas  = new  std::vector<float>;  
  m_stub_cor     = new  std::vector<float>;  
  m_stub_tp      = new  std::vector<int>;  
  m_stub_pdg     = new  std::vector<int>;  
  m_stub_pid     = new  std::vector<int>;  
  m_stub_rank    = new  std::vector<int>;  

  StubExtractor::reset();

  if (doTree)
  {
    m_OK = true;

    m_tree      = new TTree("TkStubs","Official stub info") ;

    // Branches definition

    m_tree->Branch("L1Tkevt", &n_tot_evt); // Simple evt number or event ID

    // If we don't request only matched stubs, we keep all the info
    // otherwise we skim the data file (useful for BANK generation)

    m_tree->Branch("L1TkCLUS_n",         &m_clus);
    m_tree->Branch("L1TkCLUS_x",         &m_clus_x);
    m_tree->Branch("L1TkCLUS_y",         &m_clus_y);
    m_tree->Branch("L1TkCLUS_z",         &m_clus_z);
    m_tree->Branch("L1TkCLUS_charge",    &m_clus_e);
    m_tree->Branch("L1TkCLUS_layer",     &m_clus_layer);
    m_tree->Branch("L1TkCLUS_module",    &m_clus_module);
    m_tree->Branch("L1TkCLUS_ladder",    &m_clus_ladder);
    m_tree->Branch("L1TkCLUS_seg",       &m_clus_seg);
    m_tree->Branch("L1TkCLUS_strip",     &m_clus_strip);
    m_tree->Branch("L1TkCLUS_nstrip",    &m_clus_nstrips);
    m_tree->Branch("L1TkCLUS_nsat",      &m_clus_sat);
    m_tree->Branch("L1TkCLUS_match",     &m_clus_matched);
    m_tree->Branch("L1TkCLUS_PS",        &m_clus_PS);
    m_tree->Branch("L1TkCLUS_nrows",     &m_clus_nrows);
    m_tree->Branch("L1TkCLUS_pdgID",     &m_clus_pdgID);
    m_tree->Branch("L1TkCLUS_ptGEN",     &m_clus_ptGEN);
    m_tree->Branch("L1TkCLUS_process",   &m_clus_pid);
    m_tree->Branch("L1TkCLUS_rank_PIX",  &m_clus_rank_PIX);
    m_tree->Branch("L1TkCLUS_rank_STR",  &m_clus_rank_STR);
    m_tree->Branch("L1TkCLUS_rank_PIX_4",  &m_clus_rank_PIX_4);
    m_tree->Branch("L1TkCLUS_rank_STR_4",  &m_clus_rank_STR_4);
    m_tree->Branch("L1TkCLUS_rank_PIX_8",  &m_clus_rank_PIX_8);
    m_tree->Branch("L1TkCLUS_rank_STR_8",  &m_clus_rank_STR_8);
    m_tree->Branch("L1TkCLUS_PIX",       &m_clus_pix);
    m_tree->Branch("L1TkCLUS_MULT",       &m_clus_mult);
    m_tree->Branch("L1TkCLUS_tp",        &m_clus_tp);

    m_tree->Branch("L1TkSTUB_clust1",    &m_stub_clust1);
    m_tree->Branch("L1TkSTUB_clust2",    &m_stub_clust2);
    m_tree->Branch("L1TkSTUB_cor",       &m_stub_cor);
    m_tree->Branch("L1TkSTUB_PHI0",      &m_stub_PHI0);
    m_tree->Branch("L1TkSTUB_tp",        &m_stub_tp);
    m_tree->Branch("L1TkSTUB_pdgID",     &m_stub_pdg);
    m_tree->Branch("L1TkSTUB_process",   &m_stub_pid);
    m_tree->Branch("L1TkSTUB_n",         &m_stub);
    m_tree->Branch("L1TkSTUB_pt",        &m_stub_pt);
    m_tree->Branch("L1TkSTUB_pxGEN",     &m_stub_pxGEN);
    m_tree->Branch("L1TkSTUB_pyGEN",     &m_stub_pyGEN);
    m_tree->Branch("L1TkSTUB_etaGEN",    &m_stub_etaGEN);
    m_tree->Branch("L1TkSTUB_layer",     &m_stub_layer);
    m_tree->Branch("L1TkSTUB_module",    &m_stub_module);
    m_tree->Branch("L1TkSTUB_ladder",    &m_stub_ladder);
    m_tree->Branch("L1TkSTUB_seg",       &m_stub_seg);
    m_tree->Branch("L1TkSTUB_chip",      &m_stub_chip);
    m_tree->Branch("L1TkSTUB_strip",     &m_stub_strip);
    m_tree->Branch("L1TkSTUB_x",         &m_stub_x);
    m_tree->Branch("L1TkSTUB_y",         &m_stub_y);
    m_tree->Branch("L1TkSTUB_z",         &m_stub_z);
    m_tree->Branch("L1TkSTUB_deltas",    &m_stub_deltas);
    m_tree->Branch("L1TkSTUB_X0",        &m_stub_X0);
    m_tree->Branch("L1TkSTUB_Y0",        &m_stub_Y0);
    m_tree->Branch("L1TkSTUB_Z0",        &m_stub_Z0);
    m_tree->Branch("L1TkSTUB_rank",      &m_stub_rank);
  }
}

StubExtractor::StubExtractor(TFile *a_file)
{
  std::cout << "StubExtractor object is retrieved" << std::endl;
 
  m_clus_x       = new  std::vector<float>;
  m_clus_y       = new  std::vector<float>;
  m_clus_z       = new  std::vector<float>; 
  m_clus_e       = new  std::vector<float>;
  m_clus_layer   = new  std::vector<int>;
  m_clus_module  = new  std::vector<int>;  
  m_clus_ladder  = new  std::vector<int>;  
  m_clus_seg     = new  std::vector<int>; 
  m_clus_strip   = new  std::vector<float>; 
  m_clus_used    = new  std::vector<int>;  
  m_clus_sat     = new  std::vector<int>; 
  m_clus_nstrips = new  std::vector<int>; 
  m_clus_matched = new  std::vector<int>;  
  m_clus_PS      = new  std::vector<int>;  
  m_clus_nrows   = new  std::vector<int>;
  m_clus_pid     = new  std::vector<int>;  
  m_clus_pdgID   = new  std::vector<int>;  
  m_clus_ptGEN   = new  std::vector<float>; 
  m_clus_rank_PIX= new  std::vector<int>; 
  m_clus_rank_STR= new  std::vector<int>; 
  m_clus_rank_PIX_4= new  std::vector<int>; 
  m_clus_rank_STR_4= new  std::vector<int>; 
  m_clus_rank_PIX_8= new  std::vector<int>; 
  m_clus_rank_STR_8= new  std::vector<int>; 
  m_clus_tp      = new  std::vector<int>;  
  m_clus_pix     = new  std::vector<std::vector<int> >; 
  m_clus_mult    = new  std::vector<std::vector<int> >; 

  m_stub_pt      = new  std::vector<float>;
  m_stub_ptMC    = new  std::vector<float>; 
  m_stub_pxGEN   = new  std::vector<float>; 
  m_stub_pyGEN   = new  std::vector<float>; 
  m_stub_etaGEN  = new  std::vector<float>; 
  m_stub_X0      = new  std::vector<float>; 
  m_stub_Y0      = new  std::vector<float>; 
  m_stub_Z0      = new  std::vector<float>; 
  m_stub_PHI0    = new  std::vector<float>; 
  m_stub_layer   = new  std::vector<int>; 
  m_stub_module  = new  std::vector<int>;  
  m_stub_ladder  = new  std::vector<int>; 
  m_stub_seg     = new  std::vector<int>;  
  m_stub_strip   = new  std::vector<float>; 
  m_stub_chip    = new  std::vector<int>; 
  m_stub_x       = new  std::vector<float>;  
  m_stub_y       = new  std::vector<float>;  
  m_stub_z       = new  std::vector<float>;  
  m_stub_clust1  = new  std::vector<int>;  
  m_stub_clust2  = new  std::vector<int>;  
  m_stub_cw1     = new  std::vector<int>;  
  m_stub_cw2     = new  std::vector<int>;  
  m_stub_deltas  = new  std::vector<float>;  
  m_stub_cor     = new  std::vector<float>;  
  m_stub_tp      = new  std::vector<int>;  
  m_stub_pdg     = new  std::vector<int>;  
  m_stub_pid     = new  std::vector<int>; 
  m_stub_rank    = new  std::vector<int>; 

  // Tree definition
  m_OK = false;


  StubExtractor::reset();

  m_tree = dynamic_cast<TTree*>(a_file->Get("TkStubs"));

  if (!m_tree)
  {
    std::cout << "This tree (TkStubs) doesn't exist!!!" << std::endl;
    return;
  }

  m_OK = true;
  m_matching = false;

  m_n_events = m_tree->GetEntries();

  m_tree->SetBranchAddress("L1Tkevt", &n_tot_evt); // Simple evt number or event ID
  
  // If we don't request only matched stubs, we keep all the info
  // otherwise we skim the data file (useful for BANK generation)
  
  m_tree->SetBranchAddress("L1TkCLUS_n",         &m_clus);
  m_tree->SetBranchAddress("L1TkCLUS_x",         &m_clus_x);
  m_tree->SetBranchAddress("L1TkCLUS_y",         &m_clus_y);
  m_tree->SetBranchAddress("L1TkCLUS_z",         &m_clus_z);
  m_tree->SetBranchAddress("L1TkCLUS_charge",    &m_clus_e);
  m_tree->SetBranchAddress("L1TkCLUS_layer",     &m_clus_layer);
  m_tree->SetBranchAddress("L1TkCLUS_module",    &m_clus_module);
  m_tree->SetBranchAddress("L1TkCLUS_ladder",    &m_clus_ladder);
  m_tree->SetBranchAddress("L1TkCLUS_seg",       &m_clus_seg);
  m_tree->SetBranchAddress("L1TkCLUS_strip",     &m_clus_strip);
  m_tree->SetBranchAddress("L1TkCLUS_nstrip",    &m_clus_nstrips);
  m_tree->SetBranchAddress("L1TkCLUS_nsat",      &m_clus_sat);
  m_tree->SetBranchAddress("L1TkCLUS_match",     &m_clus_matched);
  m_tree->SetBranchAddress("L1TkCLUS_PS",        &m_clus_PS);
  m_tree->SetBranchAddress("L1TkCLUS_nrows",     &m_clus_nrows);
  m_tree->SetBranchAddress("L1TkCLUS_pdgID",     &m_clus_pdgID);
  m_tree->SetBranchAddress("L1TkCLUS_ptGEN",     &m_clus_ptGEN);
  m_tree->SetBranchAddress("L1TkCLUS_process",   &m_clus_pid);
  m_tree->SetBranchAddress("L1TkCLUS_rank_PIX",  &m_clus_rank_PIX);
  m_tree->SetBranchAddress("L1TkCLUS_rank_STR",  &m_clus_rank_STR);
  m_tree->SetBranchAddress("L1TkCLUS_rank_PIX_4",  &m_clus_rank_PIX_4);
  m_tree->SetBranchAddress("L1TkCLUS_rank_STR_4",  &m_clus_rank_STR_4);
  m_tree->SetBranchAddress("L1TkCLUS_rank_PIX_8",  &m_clus_rank_PIX_8);
  m_tree->SetBranchAddress("L1TkCLUS_rank_STR_8",  &m_clus_rank_STR_8);
  m_tree->SetBranchAddress("L1TkCLUS_PIX",       &m_clus_pix);
  m_tree->SetBranchAddress("L1TkCLUS_MULT",       &m_clus_mult);
  m_tree->SetBranchAddress("L1TkCLUS_tp",        &m_clus_tp);

  m_tree->SetBranchAddress("L1TkSTUB_clust1",    &m_stub_clust1);
  m_tree->SetBranchAddress("L1TkSTUB_clust2",    &m_stub_clust2);
  m_tree->SetBranchAddress("L1TkSTUB_cor",       &m_stub_cor);
  m_tree->SetBranchAddress("L1TkSTUB_PHI0",      &m_stub_PHI0);
  m_tree->SetBranchAddress("L1TkSTUB_tp",        &m_stub_tp);
  m_tree->SetBranchAddress("L1TkSTUB_pdgID",     &m_stub_pdg);
  m_tree->SetBranchAddress("L1TkSTUB_process",   &m_stub_pid);
  m_tree->SetBranchAddress("L1TkSTUB_n",         &m_stub);
  m_tree->SetBranchAddress("L1TkSTUB_pt",        &m_stub_pt);
  m_tree->SetBranchAddress("L1TkSTUB_pxGEN",     &m_stub_pxGEN);
  m_tree->SetBranchAddress("L1TkSTUB_pyGEN",     &m_stub_pyGEN);
  m_tree->SetBranchAddress("L1TkSTUB_etaGEN",    &m_stub_etaGEN);
  m_tree->SetBranchAddress("L1TkSTUB_layer",     &m_stub_layer);
  m_tree->SetBranchAddress("L1TkSTUB_module",    &m_stub_module);
  m_tree->SetBranchAddress("L1TkSTUB_ladder",    &m_stub_ladder);
  m_tree->SetBranchAddress("L1TkSTUB_seg",       &m_stub_seg);
  m_tree->SetBranchAddress("L1TkSTUB_strip",     &m_stub_strip);
  m_tree->SetBranchAddress("L1TkSTUB_chip",      &m_stub_chip);
  m_tree->SetBranchAddress("L1TkSTUB_x",         &m_stub_x);
  m_tree->SetBranchAddress("L1TkSTUB_y",         &m_stub_y);
  m_tree->SetBranchAddress("L1TkSTUB_z",         &m_stub_z);
  m_tree->SetBranchAddress("L1TkSTUB_deltas",    &m_stub_deltas);
  m_tree->SetBranchAddress("L1TkSTUB_X0",        &m_stub_X0);
  m_tree->SetBranchAddress("L1TkSTUB_Y0",        &m_stub_Y0);
  m_tree->SetBranchAddress("L1TkSTUB_Z0",        &m_stub_Z0);
  m_tree->SetBranchAddress("L1TkSTUB_rank",      &m_stub_rank);

  std::cout << "This file contains " << m_n_events << " events..." << std::endl;

}



StubExtractor::~StubExtractor()
{}


void StubExtractor::init(const edm::EventSetup *setup)
{
  setup->get<TrackerDigiGeometryRecord>().get(theTrackerGeometry);
  setup->get<StackedTrackerGeometryRecord>().get(theStackedTrackerGeometry);
  theStackedGeometry = theStackedTrackerGeometry.product(); 

  /// Magnetic Field
  edm::ESHandle< MagneticField > magneticFieldHandle;
  setup->get< IdealMagneticFieldRecord >().get(magneticFieldHandle);
  const MagneticField* theMagneticField = magneticFieldHandle.product();
  mMagneticFieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();

}

//
// Method filling the main particle tree
//

void StubExtractor::writeInfo(const edm::Event *event, MCExtractor *mc) 
{
  StubExtractor::reset();
  ++n_tot_evt;    

  mc->clearTP(0.001,10000000.0);
  /// Sim Tracks and Vtx

  event->getByLabel( "g4SimHits", SimTrackHandle );
  event->getByLabel( "g4SimHits", SimVtxHandle );

  /// Track Trigger
  event->getByLabel( m_CLUS_tag, m_CLUS_name,    PixelDigiL1TkClusterHandle );
  event->getByLabel( m_STUB_tag, m_STUB_name,    PixelDigiL1TkStubHandle );

  /// Track Trigger MC Truth
  event->getByLabel( "TTClusterAssociatorFromPixelDigis", m_CLUS_name, MCTruthTTClusterHandle );
  event->getByLabel( "TTStubAssociatorFromPixelDigis", m_STUB_name,    MCTruthTTStubHandle );

  /// TrackingParticles
  event->getByLabel( "mix", "MergedTrackTruth", TrackingParticleHandle );
  event->getByLabel( "mix", "MergedTrackTruth", TrackingVertexHandle );

  int layer  = 0;
  int ladder = 0;
  int module = 0;
  int segs   = 0;
  int rows   = 0;

  std::vector<int> clus_coords;

  std::vector<int> clus_rows;
  std::vector<int> clus_cols;

  /// Go on only if there are L1TkCluster from PixelDigis
  if ( PixelDigiL1TkClusterHandle->size() > 0 )
  {
    /// Loop over L1TkClusters
    typename edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > >::const_iterator inputIter;
    typename edmNew::DetSet< TTCluster< Ref_PixelDigi_ > >::const_iterator contentIter;

    for ( inputIter = PixelDigiL1TkClusterHandle->begin();
	  inputIter != PixelDigiL1TkClusterHandle->end();
	  ++inputIter )
    {
      for ( contentIter = inputIter->begin();
	    contentIter != inputIter->end();
	    ++contentIter )
      {
	edm::Ref< edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > >, TTCluster< Ref_PixelDigi_ > > tempCluRef = 
	  edmNew::makeRefTo( PixelDigiL1TkClusterHandle, contentIter );

	bool genuineClu = MCTruthTTClusterHandle->isGenuine( tempCluRef );

	StackedTrackerDetId detIdClu( tempCluRef->getDetId() );
	GlobalPoint posClu      = theStackedGeometry->findAverageGlobalPosition( &(*tempCluRef) );
	MeasurementPoint coords = tempCluRef->findAverageLocalCoordinates();
	int    stack            = tempCluRef->getStackMember();

	clus_rows = tempCluRef->getRows();
	clus_cols = tempCluRef->getCols();
	clus_coords.clear();

	for (unsigned int i=0;i<clus_rows.size();++i) 
	{
	  clus_coords.push_back(clus_rows.at(i));
	  clus_coords.push_back(clus_cols.at(i));
	}

	++m_clus;

	m_clus_pix->push_back(clus_coords);

	m_clus_x->push_back(posClu.x());
	m_clus_y->push_back(posClu.y());
	m_clus_z->push_back(posClu.z());
	
	m_clus_seg->push_back(coords.y());
	m_clus_strip->push_back(coords.x());
	m_clus_nstrips->push_back(tempCluRef->findWidth());
	
	segs=2;
	rows=1024;

	if ( detIdClu.isBarrel() )
	{
	  layer  = detIdClu.iLayer()+4;
	  ladder = detIdClu.iPhi();
	  module = detIdClu.iZ()*2-1+stack;
	  if (layer<8 && stack==0) segs=32;	
	  if (layer<8) rows=960;	
	}
	else if ( detIdClu.isEndcap() )
	{	
	  layer  = 10+detIdClu.iZ()+abs(detIdClu.iSide()-2)*7;
	  ladder = detIdClu.iRing();
	  module = detIdClu.iPhi()*2-1+stack;
	  if (ladder<9 && stack==0) segs=32;
	  if (ladder<9) rows=960;
	}
	
	m_clus_layer->push_back(layer);
	m_clus_ladder->push_back(ladder);
	m_clus_module->push_back(module);
	m_clus_PS->push_back(segs);
	
	if (genuineClu)
	{
	  edm::Ptr< TrackingParticle > tpPtr = MCTruthTTClusterHandle->findTrackingParticlePtr( tempCluRef );

	  m_clus_matched->push_back(1);
	  m_clus_ptGEN->push_back(tpPtr->p4().pt()); 
	  m_clus_pdgID->push_back(tpPtr->pdgId());
	  m_clus_tp->push_back(mc->getMatchingTP(tpPtr->vertex().x(),tpPtr->vertex().y(),tpPtr->vertex().z(),
						 tpPtr->p4().px(),tpPtr->p4().py(),tpPtr->p4().pz()));
	}
	else
	{
	  m_clus_matched->push_back(0);
	  m_clus_ptGEN->push_back(0); 
	  m_clus_pdgID->push_back(0);
	  m_clus_tp->push_back(-1);
	}

	m_clus_sat->push_back(0);
	m_clus_nrows->push_back(rows);
	
      } /// End of Loop over L1TkClusters in the detId
    }
  } /// End of if ( PixelDigiL1TkClusterHandle->size() > 0 )
  

  // Clusters are built, now look at the occupancies

  int chip_idx;
  int chip_basis;

  std::vector<std::vector <int> > chip_mult;
  std::vector<std::vector<int> > chip_ranks;

  bool found;
  std::vector<int> chip_st;

  std::vector<int> strip_coord;

  int s_min;
  int s_max;

  int c_min;
  int c_max;

  int nums;

  int n_4;
  int n_8;

  chip_ranks.clear();
  chip_mult.clear();
  m_clus_rank_PIX->clear();
  m_clus_rank_STR->clear();

  for (unsigned int i=0;i<m_clus_layer->size();++i) 
  {
    strip_coord = m_clus_pix->at(i);

    s_min=2000;
    s_max=0;

    for (unsigned int j=0;j<strip_coord.size()/2;++j) 
    {
      if (strip_coord.at(2*j)<s_min) s_min = strip_coord.at(2*j);
      if (strip_coord.at(2*j)>s_max) s_max = strip_coord.at(2*j);
    }

    //    std::cout << "Cluster spans over " << s_min << "/" << s_max << std::endl;

    c_min = int(s_min/(m_clus_nrows->at(i)/8));
    c_max = int(s_max/(m_clus_nrows->at(i)/8));

    nums = m_clus_nrows->at(i)/8;

    //    std::cout << "This corresponds to chips " << c_min << " to " << c_max << std::endl;

    chip_basis=0;
    chip_basis+=(m_clus_layer->at(i))*1000000;
    chip_basis+=(m_clus_ladder->at(i))*10000;
    chip_basis+=int((m_clus_module->at(i)-1)/2)*100;
    chip_basis+=8*int(2*m_clus_seg->at(i)/m_clus_PS->at(i));

    for (int j=c_min;j<=c_max;++j) // Loop over chips
    {
      chip_idx=chip_basis+j;

      // Is this idx already registered?

      if (c_min==c_max)
      {
	n_4 = int((m_clus_nstrips->at(i)-1)/4+1);
	n_8 = int((m_clus_nstrips->at(i)-1)/8+1);
      }
      else if (j==c_min)
      {
	n_4 = int(((c_min+1)*nums-1-s_min)/4+1);
	n_8 = int(((c_min+1)*nums-1-s_min)/8+1);	
      }
      else if (j==c_max)
      {
	n_4 = int((s_max-c_max*nums-1)/4+1);
	n_8 = int((s_max-c_max*nums-1)/8+1);	
      }
      else if (j!=c_max && j!=c_min)
      {
	n_4 = int((nums-1)/4+1);
	n_8 = int((nums-1)/8+1);	
      }    

      found=false;

      for (unsigned int k=0;k<chip_mult.size();++k) 
      {  
	if (found) continue;

	if (chip_idx==chip_mult.at(k).at(0))
	{
	  found=true;
	  if (m_clus_PS->at(i)==32)
	  {
	    chip_mult.at(k).at(1)+=1;
	    chip_mult.at(k).at(3)+=n_4;
	    chip_mult.at(k).at(5)+=n_8;
	  }
	  else
	  {
	    chip_mult.at(k).at(2)+=1;
	    chip_mult.at(k).at(4)+=n_4;
	    chip_mult.at(k).at(6)+=n_8;
	  }
	}
      }

      if (!found)
      {
	chip_st.clear(); 
	chip_st.push_back(chip_idx); 
	chip_st.push_back(0);
	chip_st.push_back(0);
	chip_st.push_back(0);
	chip_st.push_back(0);
	chip_st.push_back(0);
	chip_st.push_back(0);
	
	if (m_clus_PS->at(i)==32)
	{
	  chip_st.at(1)+=1;
	  chip_st.at(3)+=n_4;
	  chip_st.at(5)+=n_8;
	}
	else
	{
	  chip_st.at(2)+=1;
	  chip_st.at(4)+=n_4;
	  chip_st.at(6)+=n_8;
	}

	chip_mult.push_back(chip_st);
      }
    }
  }

  for (unsigned int k=0;k<chip_mult.size();++k) 
  {  
    m_clus_mult->push_back(chip_mult.at(k));
  }

  //  std::cout << "This event has touched " << chip_mult.size() << " chips" << std::endl;

  for (unsigned int i=0;i<m_clus_layer->size();++i) 
  {
    chip_idx=0;
    chip_idx+=(m_clus_layer->at(i))*1000000;
    chip_idx+=(m_clus_ladder->at(i))*10000;
    chip_idx+=int((m_clus_module->at(i)-1)/2)*100;
    chip_idx+=8*int(2*m_clus_seg->at(i)/m_clus_PS->at(i));
    chip_idx+=int(m_clus_strip->at(i)/(m_clus_nrows->at(i)/8));

    found=false;

    for (unsigned int k=0;k<chip_mult.size();++k) 
    {  
      if (found) continue;

      if (chip_idx==chip_mult.at(k).at(0))
      {
	found=true;
	m_clus_rank_PIX->push_back(chip_mult.at(k).at(1));
	m_clus_rank_STR->push_back(chip_mult.at(k).at(2));
	m_clus_rank_PIX_4->push_back(chip_mult.at(k).at(3));
	m_clus_rank_STR_4->push_back(chip_mult.at(k).at(4));
	m_clus_rank_PIX_8->push_back(chip_mult.at(k).at(5));
	m_clus_rank_STR_8->push_back(chip_mult.at(k).at(6));
      }
    }
  }



  int clust1     = -1;
  int clust2     = -1;
 
  /// Go on only if there are L1TkStub from PixelDigis
  if ( PixelDigiL1TkStubHandle->size() > 0 )
  {
    /// Loop over L1TkStubs

    typename edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >::const_iterator otherInputIter;
    typename edmNew::DetSet< TTStub< Ref_PixelDigi_ > >::const_iterator otherContentIter;

    for ( otherInputIter = PixelDigiL1TkStubHandle->begin();
	  otherInputIter != PixelDigiL1TkStubHandle->end();
	  ++otherInputIter )
    {
      for ( otherContentIter = otherInputIter->begin();
	    otherContentIter != otherInputIter->end();
	    ++otherContentIter )
      {
	/// Make the reference to be put in the map
	edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > tempStubPtr = 
	  edmNew::makeRefTo( PixelDigiL1TkStubHandle, otherContentIter );

	StackedTrackerDetId detIdStub( tempStubPtr->getDetId() );
      
	// The stub position is the inner cluster position
	//

	GlobalPoint posStub = theStackedGeometry->findGlobalPosition( &(*tempStubPtr) );

	++m_stub;

	double displStub    = tempStubPtr->getTriggerDisplacement();
	double offsetStub   = tempStubPtr->getTriggerOffset();
	
	bool genuineStub    = MCTruthTTStubHandle->isGenuine( tempStubPtr );
    
	const GeomDetUnit* det0 = theStackedGeometry->idToDetUnit( detIdStub, 0 );

	/// Find pixel pitch and topology related information
	const PixelGeomDetUnit* pix0 = dynamic_cast< const PixelGeomDetUnit* >( det0 );
	const PixelTopology* top0    = dynamic_cast< const PixelTopology* >( &(pix0->specificTopology()) );
	
	segs = top0->ncolumns();
	rows = top0->nrows();
	
	clust1 = StubExtractor::getClust1Idx(posStub.x(),posStub.y(),posStub.z());
	clust2 = StubExtractor::getClust2Idx(clust1,m_clus_strip->at(clust1)+displStub);
	
	m_stub_x->push_back(posStub.x());
	m_stub_y->push_back(posStub.y());
	m_stub_z->push_back(posStub.z());
	m_stub_clust1->push_back(clust1);
	m_stub_clust2->push_back(clust2);
	m_stub_seg->push_back(m_clus_seg->at(clust1));
	m_stub_chip->push_back(m_clus_strip->at(clust1)/(rows/8));
	m_stub_strip->push_back(m_clus_strip->at(clust1));
	m_stub_deltas->push_back(displStub-offsetStub);
	m_stub_cor->push_back(offsetStub);
	m_stub_pt->push_back(theStackedGeometry->findRoughPt( mMagneticFieldStrength, &(*tempStubPtr) ));

	
	if ( detIdStub.isBarrel() )
	{
	  layer  = detIdStub.iLayer()+4;
	  ladder = detIdStub.iPhi()-1;
	  module = detIdStub.iZ()-1;
	}
	else if ( detIdStub.isEndcap() )
	{	
	  layer  = 10+detIdStub.iZ()+abs(detIdStub.iSide()-2)*7;
	  ladder = detIdStub.iRing()-1;
	  module = detIdStub.iPhi()-1;
	}

	m_stub_layer->push_back(layer);
	m_stub_ladder->push_back(ladder);
	m_stub_module->push_back(module);
	m_stub_rank->push_back(0);

	if ( genuineStub )
	{
	  
	  edm::Ptr< TrackingParticle > tpPtr = MCTruthTTStubHandle->findTrackingParticlePtr( tempStubPtr );
	  
	  m_stub_pxGEN->push_back(tpPtr->p4().px());
	  m_stub_pyGEN->push_back(tpPtr->p4().py());
	  m_stub_etaGEN->push_back(tpPtr->momentum().eta());
	  m_stub_pdg->push_back(tpPtr->pdgId());
	  m_stub_pid->push_back(0);
	  m_stub_X0->push_back(tpPtr->vertex().x());
	  m_stub_Y0->push_back(tpPtr->vertex().y());
	  m_stub_Z0->push_back(tpPtr->vertex().z());
	  m_stub_PHI0->push_back(tpPtr->momentum().phi());
	  m_stub_tp->push_back(mc->getMatchingTP(tpPtr->vertex().x(),tpPtr->vertex().y(),tpPtr->vertex().z(),
						 tpPtr->p4().px(),tpPtr->p4().py(),tpPtr->p4().pz()));
	}
	else
	{
	  m_stub_pxGEN->push_back(0);
	  m_stub_pyGEN->push_back(0);
	  m_stub_etaGEN->push_back(0);
	  m_stub_X0->push_back(0);
	  m_stub_Y0->push_back(0);
	  m_stub_Z0->push_back(0);
	  m_stub_PHI0->push_back(0);
	  m_stub_pdg->push_back(0);
	  m_stub_pid->push_back(0);
	  m_stub_tp->push_back(-1);
	}

	
	chip_idx=0;
	chip_idx+=(layer)*1000000;
	chip_idx+=(ladder)*10000;
	chip_idx+=(module)*100;
	chip_idx+=8*int(2*m_clus_seg->at(clust1)/m_clus_PS->at(clust1));
	chip_idx+=int(m_clus_strip->at(clust1)/(m_clus_nrows->at(clust1)/8));
	
	// Is this idx already registered?
	
	found=false;
	
	for (unsigned int k=0;k<chip_ranks.size();++k) 
	{  
	  if (found) continue;

	  if (chip_idx==chip_ranks.at(k).at(0))
	  {
	    found=true;
	    chip_ranks.at(k).push_back(m_stub_pt->size()-1);
	  }
	}

	if (!found)
	{
	  chip_st.clear(); 
	  chip_st.push_back(chip_idx); 
	  chip_st.push_back(m_stub_pt->size()-1); 
	  chip_ranks.push_back(chip_st);
	}
      } /// End of loop over L1TkStubs
    } 
  } /// End of if ( PixelDigiL1TkStubHandle->size() > 0 ) 
 
  
  std::vector<int> rank;
  std::vector<float> bends;
  std::vector<float> strip;

  float b_min;
  float s_min2;
  unsigned int   i_min;
  int tmp;

  for (unsigned int k=0;k<chip_ranks.size();++k) 
  { 
    chip_st=chip_ranks.at(k); 

    if (chip_st.size()==2) continue; // Only one stub, rank is 0

    rank.clear();
    bends.clear();
    strip.clear();

    for (unsigned int l=1;l<chip_st.size();++l) 
    { 
      rank.push_back(chip_st.at(l));
      bends.push_back(fabs(m_stub_deltas->at(chip_st.at(l))));
      strip.push_back(m_stub_strip->at(chip_st.at(l)));
    }

    // We now have the list, we do the sorting

    for (unsigned int l=0;l<rank.size()-1;++l) 
    { 
      i_min=l;
      b_min=bends.at(l);
      s_min2=strip.at(l);

      for (unsigned int ll=l+1;ll<rank.size();++ll) 
      { 
	if (bends.at(ll)<b_min)
	{
	  i_min=ll;
	  b_min=bends.at(ll);
	  s_min2=strip.at(ll);
	}
	if (bends.at(ll)==b_min && strip.at(ll)<s_min2)
	{
	  i_min=ll;
	  b_min=bends.at(ll);
	  s_min2=strip.at(ll);
	}
      }

      if (i_min!=l)
      {
	tmp=rank.at(i_min);
	rank.at(i_min)=rank.at(l);
	bends.at(i_min)=bends.at(l);
	strip.at(i_min)=strip.at(l);
	
	rank.at(l)=tmp;
	bends.at(l)=b_min;
	strip.at(l)=s_min2;
      }
    }

    for (unsigned int l=0;l<rank.size();++l) m_stub_rank->at(rank.at(l)) = l; 
  }

  StubExtractor::fillTree();
}


//
// Method getting the info from an input file
//

void StubExtractor::getInfo(int ievt) 
{
  m_tree->GetEntry(ievt); 
}

// Method initializing everything (to do for each event)

void StubExtractor::reset()
{
  m_clus = 0;
  m_stub = 0;

  m_clus_x->clear(); 
  m_clus_y->clear(); 
  m_clus_z->clear(); 
  m_clus_e->clear(); 
  m_clus_layer->clear(); 
  m_clus_module->clear();
  m_clus_ladder->clear();
  m_clus_seg->clear();   
  m_clus_strip->clear(); 
  m_clus_sat->clear();   
  m_clus_nstrips->clear();
  m_clus_used->clear();   
  m_clus_matched->clear();
  m_clus_PS->clear();
  m_clus_nrows->clear();
  m_clus_pid->clear();
  m_clus_pdgID->clear();
  m_clus_ptGEN->clear();
  m_clus_rank_PIX->clear(); 
  m_clus_rank_STR->clear(); 
  m_clus_rank_PIX_4->clear(); 
  m_clus_rank_STR_4->clear();
  m_clus_rank_PIX_8->clear(); 
  m_clus_rank_STR_8->clear();
  m_clus_pix->clear(); 
  m_clus_tp->clear(); 
  m_clus_mult->clear(); 

  m_stub_X0->clear();     
  m_stub_Y0->clear();     
  m_stub_Z0->clear();     
  m_stub_PHI0->clear();     
  m_stub_tp->clear();     
  m_stub_pt->clear();     
  m_stub_ptMC->clear();   
  m_stub_pxGEN->clear();  
  m_stub_pyGEN->clear();  
  m_stub_etaGEN->clear();  
  m_stub_layer->clear();  
  m_stub_module->clear(); 
  m_stub_ladder->clear(); 
  m_stub_seg->clear();  
  m_stub_chip->clear();   
  m_stub_strip->clear(); 
  m_stub_x->clear(); 
  m_stub_y->clear(); 
  m_stub_z->clear(); 
  m_stub_clust1->clear(); 
  m_stub_clust2->clear(); 
  m_stub_cw1->clear(); 
  m_stub_cw2->clear(); 
  m_stub_deltas->clear(); 
  m_stub_cor->clear(); 
  m_stub_pdg->clear();
  m_stub_pid->clear();
  m_stub_rank->clear(); 

}


void StubExtractor::fillTree()
{
  m_tree->Fill(); 
}
 
void StubExtractor::fillSize(int size)
{
  m_clus=size;
}

int  StubExtractor::getSize()
{
  return m_clus;
}

int  StubExtractor::getClust1Idx(float x, float y, float z)
{
  for (int i=0;i<m_clus;++i) // Loop over clusters
  { 
    if (fabs(m_clus_x->at(i)-x) > 0.001) continue;
    if (fabs(m_clus_y->at(i)-y) > 0.001) continue;
    if (fabs(m_clus_z->at(i)-z) > 0.001) continue;

    return i;
  }

  std::cout << "STRANGE: a stub without matching 1 cluster..." << std::endl;

  return -1;
}

int  StubExtractor::getClust2Idx(int idx1, float dist)
{
  float dmax = 20;
  int idx2 = -1;

  for (int i=0;i<m_clus;++i) // Loop over clusters
  { 
    if (m_clus_layer->at(i)!=m_clus_layer->at(idx1)) continue;
    if (m_clus_ladder->at(i)!=m_clus_ladder->at(idx1)) continue;
    if (m_clus_module->at(i)-1!=m_clus_module->at(idx1)) continue;

    //    std::cout  << m_clus_layer->at(idx1) << " / " <<  m_clus_ladder->at(idx1) << " / " 
    //               << m_clus_module->at(idx1) << " / " <<  m_clus_module->at(i) << " / " 
    //               << m_clus_strip->at(idx1) << " / " <<  m_clus_strip->at(i) << " / " 
    //    	       << fabs(m_clus_strip->at(i)-dist) << " / " << dist << std::endl;

    if (fabs(m_clus_strip->at(i)-dist)<dmax)
    {  
      dmax = fabs(m_clus_strip->at(i)-dist);
      idx2 = i;
    }
  }

  if (idx2==-1)
    std::cout << "STRANGE: a stub without matching 2 cluster..." << std::endl;

  return idx2;
}

int  StubExtractor::get_id(int lay,int lad,int mod,float x,float y,float z)
{
  int idx = -1;

  for (int i=0;i<m_stub;++i) // Loop over stubs
  { 
    if (idx!=-1) return idx;

    if (m_stub_layer->at(i)  !=lay) continue;
    if (m_stub_ladder->at(i) !=lad) continue;
    if (m_stub_module->at(i) !=mod) continue;
    if (fabs(m_stub_x->at(i)-x)>0.001) continue;
    if (fabs(m_stub_y->at(i)-y)>0.001) continue;
    if (fabs(m_stub_z->at(i)-z)>0.001) continue;

    idx=i;
  }

  if (idx!=-1) return idx;

  std::cout << "StubExtractor::get_id: you are not supposed to get here!" << std::endl;

  return idx;
}
