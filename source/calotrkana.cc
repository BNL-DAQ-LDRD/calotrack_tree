//____________________________________________________________________________..
//
//
//____________________________________________________________________________..

#include "calotrkana.h"

// Fun4all stuff
#include <ffaobjects/EventHeader.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <phool/PHCompositeNode.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

// ROOT stuff
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TTree.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

//just in case in the future I want to get HEPMC info
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h> // for GenVertex, GenVertex::part...
#pragma GCC diagnostic pop

#include <HepMC/GenParticle.h> // for GenParticle
#include <HepMC/GenRanges.h>
#include <HepMC/HeavyIon.h>      // for HeavyIon
#include <HepMC/IteratorRange.h> // for children, descendants
#include <HepMC/SimpleVector.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

//tracking stuff
#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>  // for getTrkrId, getHit...
#include <trackbase/TrkrHitTruthAssoc.h>

//____________________________________________________________________________..
calotrkana::calotrkana(const std::string &name = "calotrkana",
               const std::string &outName = "ZDCresult")
    : SubsysReco(name), Outfile(outName) {
  std::cout << "calotrkana::calotrkana(const std::string &name) Calling ctor"
            << std::endl;
}

//____________________________________________________________________________..
calotrkana::~calotrkana() {
  std::cout << "calotrkana::~calotrkana() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int calotrkana::Init(PHCompositeNode *topNode) {

  std::cout << "calotrkana::Init(PHCompositeNode *topNode) Initializing"
            << std::endl;
  out = new TFile(Outfile.c_str(), "RECREATE");
  T = new TTree("T", "Tree for Reco Hits and truth");
  T->Branch("nHits", &m_nHits, "nHits/I");
  T->Branch("Hit_E", &m_Hit_E, "Hit_E[nHits]/F");
  T->Branch("Hit_x", &m_Hit_x, "Hit_x[nHits]/F");
  T->Branch("Hit_y", &m_Hit_y, "Hit_y[nHits]/F");
  T->Branch("Hit_z", &m_Hit_z, "Hit_z[nHits]/F");
  T->Branch("Hit_t", &m_Hit_t, "Hit_t[nHits]/F");

  T->Branch("nParticles", &m_nParticles, "nParticles/I");
  T->Branch("particle_pid", &m_particle_pid, "particle_pid[nParticles]/I");
  T->Branch("particle_energy", &m_particle_energy, "particle_energy[nParticles]/F");
  T->Branch("particle_px", &m_particle_px, "particle_px[nParticles]/F");
  T->Branch("particle_py", &m_particle_py, "particle_py[nParticles]/F");
  T->Branch("particle_pz", &m_particle_pz, "particle_pz[nParticles]/F");
  T->Branch("particle_vtx_x", &m_particle_vtx_x, "particle_vtx_x[nParticles]/F");
  T->Branch("particle_vtx_y", &m_particle_vtx_y, "particle_vtx_y[nParticles]/F");
  T->Branch("particle_vtx_z", &m_particle_vtx_z, "particle_vtx_z[nParticles]/F");


  Fun4AllServer *se = Fun4AllServer::instance();
  se->Print("NODETREE");


  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int calotrkana::InitRun(PHCompositeNode *topNode) {
  std::cout
      << "calotrkana::InitRun(PHCompositeNode *topNode) Initializing for Run XXX"
      << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int calotrkana::process_event(PHCompositeNode *topNode) {

  // truth container here
  PHG4TruthInfoContainer *truthinfo =
      findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truthinfo) {
    std::cout << "calotrkana::process_event(PHCompositeNode *topNode) No truth "
                 "info found"
              << std::endl;
  }

  // loop over truth primary particles
  PHG4TruthInfoContainer::ConstRange range =
      truthinfo->GetPrimaryParticleRange();
  // set of primary pi0 trk id
  std::set<int> primary_pi0_trk_id;
  for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first;
       truth_itr != range.second; ++truth_itr) {
    PHG4Particle *truth = truth_itr->second;
    // get track id and check if embedded
    int vtxid = truth->get_vtx_id();
    PHG4VtxPoint *vtx = truthinfo->GetVtx(vtxid);

    m_particle_pid[m_nParticles] = truth->get_pid();
    m_particle_energy[m_nParticles] = truth->get_e();
    m_particle_px[m_nParticles] = truth->get_px();
    m_particle_py[m_nParticles] = truth->get_py();
    m_particle_pz[m_nParticles] = truth->get_pz();
    m_particle_vtx_x[m_nParticles] = vtx->get_x();
    m_particle_vtx_y[m_nParticles] = vtx->get_y();
    m_particle_vtx_z[m_nParticles] = vtx->get_z();
    m_nParticles++;
    if(m_nParticles >= ptruthmaxlength){
      std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nParticles exceeds max length" << std::endl;
      exit(1);
    }

    
  }

  // CEMC
  TowerInfoContainer *CEMC_towers_sim =
      findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  RawTowerGeomContainer *CEMC_geom =
      findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");

  int CEMCsize = CEMC_towers_sim->size();
  for (int i = 0; i < CEMCsize; i++) {
    TowerInfo *tower = CEMC_towers_sim->get_tower_at_channel(i);

    if (tower->get_energy() < 0.05)
      continue;
    unsigned int towerkey = CEMC_towers_sim->encode_key(i);
    int ieta = CEMC_towers_sim->getTowerEtaBin(towerkey);
    int iphi = CEMC_towers_sim->getTowerPhiBin(towerkey);
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(
        RawTowerDefs::convert_name_to_caloid("CEMC"), ieta, iphi);
    float tower_x = CEMC_geom->get_tower_geometry(key)->get_center_x();
    float tower_y = CEMC_geom->get_tower_geometry(key)->get_center_y();
    float tower_z = CEMC_geom->get_tower_geometry(key)->get_center_z();

    m_Hit_E[m_nHits] = tower->get_energy();
    m_Hit_x[m_nHits] = tower_x;
    m_Hit_y[m_nHits] = tower_y;
    m_Hit_z[m_nHits] = tower_z;
    m_Hit_t[m_nHits] = tower->get_time_float();
    m_nHits++;

    if(m_nHits >= recomaxlength){
      std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nHits exceeds max length" << std::endl;
      exit(1);
    }
  }

  // IHCAL
  TowerInfoContainer *IHCAL_towers_sim =
      findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  RawTowerGeomContainer *IHCAL_geom =
      findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");

  int IHCALsize = IHCAL_towers_sim->size();
  for (int i = 0; i < IHCALsize; i++) {
    TowerInfo *tower = IHCAL_towers_sim->get_tower_at_channel(i);

    if (tower->get_energy() < 0.05)
      continue;
    unsigned int towerkey = IHCAL_towers_sim->encode_key(i);
    int ieta = IHCAL_towers_sim->getTowerEtaBin(towerkey);
    int iphi = IHCAL_towers_sim->getTowerPhiBin(towerkey);
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(
        RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
    float tower_x = IHCAL_geom->get_tower_geometry(key)->get_center_x();
    float tower_y = IHCAL_geom->get_tower_geometry(key)->get_center_y();
    float tower_z = IHCAL_geom->get_tower_geometry(key)->get_center_z();

    m_Hit_E[m_nHits] = tower->get_energy();
    m_Hit_x[m_nHits] = tower_x;
    m_Hit_y[m_nHits] = tower_y;
    m_Hit_z[m_nHits] = tower_z;
    m_Hit_t[m_nHits] = tower->get_time_float();
    m_nHits++;

    if(m_nHits >= recomaxlength){
      std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nHits exceeds max length" << std::endl;
      exit(1);
    }
  }

  // OHCAL
  TowerInfoContainer *OHCAL_towers_sim =
      findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  RawTowerGeomContainer *OHCAL_geom =
      findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  if (!OHCAL_towers_sim || !OHCAL_geom) {
    std::cout << "calotrkana::process_event(PHCompositeNode *topNode) No OHCAL "
                 "towers found"
              << std::endl;
  } 
  int OHCALsize = OHCAL_towers_sim->size();
  for (int i = 0; i < OHCALsize; i++) {
    TowerInfo *tower = OHCAL_towers_sim->get_tower_at_channel(i);

    if (tower->get_energy() < 0.05)
      continue;
    unsigned int towerkey = OHCAL_towers_sim->encode_key(i);
    int ieta = OHCAL_towers_sim->getTowerEtaBin(towerkey);
    int iphi = OHCAL_towers_sim->getTowerPhiBin(towerkey);
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(
        RawTowerDefs::CalorimeterId::HCALOUT, ieta, iphi);
    float tower_x = OHCAL_geom->get_tower_geometry(key)->get_center_x();
    float tower_y = OHCAL_geom->get_tower_geometry(key)->get_center_y();
    float tower_z = OHCAL_geom->get_tower_geometry(key)->get_center_z();

    m_Hit_E[m_nHits] = tower->get_energy();
    m_Hit_x[m_nHits] = tower_x;
    m_Hit_y[m_nHits] = tower_y;
    m_Hit_z[m_nHits] = tower_z;
    m_Hit_t[m_nHits] = tower->get_time_float();
    m_nHits++;

    if(m_nHits >= recomaxlength){
      std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nHits exceeds max length" << std::endl;
      exit(1);
    }
  }

  // tracking
  ActsGeometry *m_tGeometry =
      findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry) {
    std::cout << PHWHERE << "No acts tracking geometry, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  TrkrClusterContainer *clustermap =
      findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!clustermap) {
    std::cout << "calotrkana::process_event(PHCompositeNode *topNode) No cluster "
                 "map found"
              << std::endl;
  }
  for (const auto &hitsetkey :
       clustermap->getHitSetKeys()) {
    auto range = clustermap->getClusters(hitsetkey);
    for (auto clusterIter = range.first; clusterIter != range.second; ++clusterIter){
      const auto &key = clusterIter->first;
      const auto &cluster = clusterIter->second;

      const auto global = m_tGeometry->getGlobalPosition(key, cluster);

      float x = global(0);
      float y = global(1);
      float z = global(2);

      float e = cluster->getAdc();
      float t = cluster->getTime();

      m_Hit_E[m_nHits] = e;
      m_Hit_x[m_nHits] = x;
      m_Hit_y[m_nHits] = y;
      m_Hit_z[m_nHits] = z;
      m_Hit_t[m_nHits] = t;
      m_nHits++;

      if(m_nHits >= recomaxlength){
        std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nHits exceeds max length" << std::endl;
        exit(1);
      }
    }
  }

  T->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int calotrkana::ResetEvent(PHCompositeNode *topNode) {

  m_nHits = 0;
  m_nParticles = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________________________________________________..
int calotrkana::End(PHCompositeNode *topNode) {
  out->cd();

  T->Write();

  out->Close();
  delete out;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int calotrkana::Reset(PHCompositeNode *topNode) {
  std::cout << "calotrkana::Reset(PHCompositeNode *topNode) being Reset"
            << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void calotrkana::Print(const std::string &what) const {
  std::cout << "calotrkana::Print(const std::string &what) const Printing info for "
            << what << std::endl;
}
