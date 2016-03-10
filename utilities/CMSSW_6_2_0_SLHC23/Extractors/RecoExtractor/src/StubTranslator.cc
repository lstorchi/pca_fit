#include "../interface/StubTranslator.h"

StubTranslator::StubTranslator()
{
  std::cout << "Entering OffStub Translator" << std::endl;

  m_stub_etaGEN  = new  std::vector<float>; 
  m_stub_strip   = new  std::vector<float>; 
  m_stub_ptGEN   = new  std::vector<float>;  
  m_stub_x       = new  std::vector<float>;  
  m_stub_y       = new  std::vector<float>;  
  m_stub_z       = new  std::vector<float>;
  m_stub_bend    = new  std::vector<float>;    
  m_stub_modid   = new  std::vector<int>;  
  m_stub_X0      = new  std::vector<float>;  
  m_stub_Y0      = new  std::vector<float>;  
  m_stub_Z0      = new  std::vector<float>;
  m_stub_PHI0    = new  std::vector<float>;
  m_stub_pdg     = new  std::vector<int>;

  StubTranslator::reset();

  // Then the rootuple

  m_tree_L1TrackTrigger = new TTree("BankStubs","Stubs for bank");  

  m_tree_L1TrackTrigger->Branch("STUB_n",           &m_stub);
  m_tree_L1TrackTrigger->Branch("STUB_ptGEN",       &m_stub_ptGEN);
  m_tree_L1TrackTrigger->Branch("STUB_etaGEN",      &m_stub_etaGEN);
  m_tree_L1TrackTrigger->Branch("STUB_x",           &m_stub_x);
  m_tree_L1TrackTrigger->Branch("STUB_y",           &m_stub_y);
  m_tree_L1TrackTrigger->Branch("STUB_z",           &m_stub_z);
  m_tree_L1TrackTrigger->Branch("STUB_bend",        &m_stub_bend);
  m_tree_L1TrackTrigger->Branch("STUB_modid",       &m_stub_modid);
  m_tree_L1TrackTrigger->Branch("STUB_strip",       &m_stub_strip);
  m_tree_L1TrackTrigger->Branch("STUB_X0",          &m_stub_X0);
  m_tree_L1TrackTrigger->Branch("STUB_Y0",          &m_stub_Y0);
  m_tree_L1TrackTrigger->Branch("STUB_Z0",          &m_stub_Z0);
  m_tree_L1TrackTrigger->Branch("STUB_PHI0",        &m_stub_PHI0);
  m_tree_L1TrackTrigger->Branch("STUB_pdg" ,        &m_stub_pdg);


}

StubTranslator::~StubTranslator()
{}

void StubTranslator::reset()
{
  m_stub = 0;

  m_stub_etaGEN->clear();  
  m_stub_strip->clear(); 
  m_stub_ptGEN->clear();
  m_stub_x->clear();
  m_stub_y->clear();
  m_stub_z->clear();
  m_stub_bend->clear();
  m_stub_modid->clear();  
  m_stub_X0->clear();
  m_stub_Y0->clear();
  m_stub_Z0->clear();
  m_stub_PHI0->clear();
  m_stub_pdg->clear();

}

void StubTranslator::do_translation(StubExtractor *st)
{

  int n_stubs = st->getNStub();

  StubTranslator::reset();

  for (int i=0;i<n_stubs;++i)
  {
    if (st->getStub_tp(i)!=0) continue;

    ++m_stub;
    
    m_stub_etaGEN->push_back(st->getStub_etaGen(i));  
    m_stub_strip->push_back(st->getStub_strip(i)); 
    m_stub_ptGEN->push_back(st->getStub_ptGen(i));
    m_stub_modid->push_back(st->getStub_ID(i)); 
    m_stub_x->push_back(st->getStub_x(i)); 
    m_stub_y->push_back(st->getStub_y(i)); 
    m_stub_z->push_back(st->getStub_z(i)); 
    m_stub_bend->push_back(st->getStub_bend(i)); 
    m_stub_X0->push_back(st->getStub_X0(i)); 
    m_stub_Y0->push_back(st->getStub_Y0(i)); 
    m_stub_Z0->push_back(st->getStub_Z0(i)); 
    m_stub_PHI0->push_back(st->getStub_PHI0(i));
    m_stub_pdg->push_back(st->getStub_pdg(i));

  }
}

// Fill the root tree containing analysis results

void StubTranslator::fillTree()
{
  m_tree_L1TrackTrigger->Fill(); 
}
