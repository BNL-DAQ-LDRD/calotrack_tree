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

#include <mbd/MbdGeom.h>
#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <epd/EpdGeom.h>

//tracking stuff
#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>  // for getTrkrId, getHit...
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/InttDefs.h>

#include <g4eval/SvtxEvalStack.h>

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
  T->Branch("Hit_detid", &m_Hit_detid, "Hit_detid[nHits]/I");
  T->Branch("nParticles", &m_nParticles, "nParticles/I");
  T->Branch("particle_pid", &m_particle_pid, "particle_pid[nParticles]/I");
  T->Branch("particle_energy", &m_particle_energy, "particle_energy[nParticles]/F");
  T->Branch("particle_px", &m_particle_px, "particle_px[nParticles]/F");
  T->Branch("particle_py", &m_particle_py, "particle_py[nParticles]/F");
  T->Branch("particle_pz", &m_particle_pz, "particle_pz[nParticles]/F");
  T->Branch("particle_vtx_x", &m_particle_vtx_x, "particle_vtx_x[nParticles]/F");
  T->Branch("particle_vtx_y", &m_particle_vtx_y, "particle_vtx_y[nParticles]/F");
  T->Branch("particle_vtx_z", &m_particle_vtx_z, "particle_vtx_z[nParticles]/F");
  T->Branch("particle_track_id", &m_particle_track_id, "particle_track_id[nParticles]/I");
  T->Branch("particle_primary_id", &m_particle_primary_id, "particle_primary_id[nParticles]/I");
  T->Branch("nRecoClusters", &m_nRecoClusters, "nRecoClusters/I");
  T->Branch("reco_cluster_E", &m_reco_cluster_E, "reco_cluster_E[nRecoClusters]/F");
  T->Branch("reco_cluster_x", &m_reco_cluster_x, "reco_cluster_x[nRecoClusters]/F");
  T->Branch("reco_cluster_y", &m_reco_cluster_y, "reco_cluster_y[nRecoClusters]/F");
  T->Branch("reco_cluster_z", &m_reco_cluster_z, "reco_cluster_z[nRecoClusters]/F");
  T->Branch("reco_cluster_t", &m_reco_cluster_t, "reco_cluster_t[nRecoClusters]/F");
  T->Branch("reco_cluster_detid", &m_reco_cluster_detid, "reco_cluster_detid[nRecoClusters]/I");
  T->Branch("reco_cluster_trcluster_id", &m_reco_cluster_trcluster_id, "reco_cluster_trcluster_id[nRecoClusters]/i");
  T->Branch("nTruthClusters", &m_nTruthClusters, "nTruthClusters/I");
  T->Branch("truth_cluster_E", &m_truth_cluster_E, "truth_cluster_E[nTruthClusters]/F");
  T->Branch("truth_cluster_x", &m_truth_cluster_x, "truth_cluster_x[nTruthClusters]/F");
  T->Branch("truth_cluster_y", &m_truth_cluster_y, "truth_cluster_y[nTruthClusters]/F");
  T->Branch("truth_cluster_z", &m_truth_cluster_z, "truth_cluster_z[nTruthClusters]/F");
  T->Branch("truth_cluster_t", &m_truth_cluster_t, "truth_cluster_t[nTruthClusters]/F");
  T->Branch("truth_cluster_detid", &m_truth_cluster_detid, "truth_cluster_detid[nTruthClusters]/I");
  T->Branch("truth_cluster_id", &m_truth_cluster_id, "truth_cluster_id[nTruthClusters]/i");
  T->Branch("truth_cluster_trparticle_track_id", &m_truth_cluster_trparticle_track_id, "truth_cluster_trparticle_track_id[nTruthClusters]/I");



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
  // set of primary particles
  std::set<PHG4Particle*> primary_particles;
  for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first;
       truth_itr != range.second; ++truth_itr) {
    PHG4Particle *truth = truth_itr->second;
    primary_particles.insert(truth);

    
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
    m_Hit_t[m_nHits] = tower->get_time_float()*sampletons;
    m_Hit_detid[m_nHits] = static_cast<int>(calotrkana::cemcId);
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
    m_Hit_t[m_nHits] = tower->get_time_float()*sampletons;
    m_Hit_detid[m_nHits] = static_cast<int>(calotrkana::ihcalId);
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
    m_Hit_t[m_nHits] = tower->get_time_float()*sampletons;
    m_Hit_detid[m_nHits] = static_cast<int>(calotrkana::ohcalId);
    m_nHits++;

    if(m_nHits >= recomaxlength){
      std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nHits exceeds max length" << std::endl;
      exit(1);
    }
  }

  //SEPD
  TowerInfoContainer *EPD_towers_sim =
      findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_EPD");
  EpdGeom *EPD_geom =
      findNode::getClass<EpdGeom>(topNode, "TOWERGEOM_EPD");
  if (!OHCAL_towers_sim || !OHCAL_geom) {
    std::cout << "calotrkana::process_event(PHCompositeNode *topNode) No EPD "
                 "towers found"
              << std::endl;
  }

  unsigned int EPDsize = EPD_towers_sim->size();
      for (unsigned int i = 0; i < EPDsize; i++)
      {
        TowerInfo *tower = EPD_towers_sim ->get_tower_at_channel(i);
        unsigned int key = TowerInfoDefs::encode_epd(i);
        if (tower->get_energy() < 0.05)
          continue;
        float r = EPD_geom->get_r(key);
        float phi = EPD_geom->get_phi(key);
        float z = EPD_geom->get_z(key);

        float x = r * cos(phi);
        float y = r * sin(phi);

        m_Hit_E[m_nHits] = tower->get_energy();
        m_Hit_x[m_nHits] = x;
        m_Hit_y[m_nHits] = y;
        m_Hit_z[m_nHits] = z;
        m_Hit_t[m_nHits] = tower->get_time();
        m_Hit_detid[m_nHits] = static_cast<int>(calotrkana::epdId);
        m_nHits++;

        if(m_nHits >= recomaxlength){
          std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nHits exceeds max length" << std::endl;
          exit(1);
        }

      }

      //MBD
      MbdPmtContainer *mbdpmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
      MbdGeom *mbdgeom = findNode::getClass<MbdGeom>(topNode, "MbdGeom");
      if (!mbdpmts || !mbdgeom) {
        std::cout << "calotrkana::process_event(PHCompositeNode *topNode) No MBD "
                     "pmts found"
                  << std::endl;
      }
      for (int ipmt = 0; ipmt < mbdpmts->get_npmt(); ipmt++)
      {
        float mbdq = mbdpmts->get_pmt(ipmt)->get_q();
        if (mbdq < 0.25) //from MBD reco hit threshold
          continue;
        float mbdt = mbdpmts->get_pmt(ipmt)->get_time();
        //check if mbdt is undefined
        if(mbdt > 1E6)
          continue;

        float x = mbdgeom->get_x(ipmt);
        float y = mbdgeom->get_y(ipmt);
        float z = mbdgeom->get_z(ipmt);

        m_Hit_E[m_nHits] = mbdq;
        m_Hit_x[m_nHits] = x;
        m_Hit_y[m_nHits] = y;
        m_Hit_z[m_nHits] = z;
        m_Hit_t[m_nHits] = mbdt;
        m_Hit_detid[m_nHits] = static_cast<int>(calotrkana::mbdId);
        m_nHits++;

        if(m_nHits >= recomaxlength){
          std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nHits exceeds max length" << std::endl;
          exit(1);
        }
      }
  //find secondary particles associated with tack clusters
  if (!m_svtxEvalStack)
  {
    m_svtxEvalStack = new SvtxEvalStack(topNode);
  }
  //do eval
  m_svtxEvalStack->next_event(topNode);

  /// Get the cluster evaluator
  SvtxClusterEval *clustereval = m_svtxEvalStack->get_cluster_eval();
  //clustereval->set_verbosity(1);
  SvtxTruthEval *trutheval = m_svtxEvalStack->get_truth_eval();
  //trutheval->set_verbosity(2);
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
  std::set<PHG4Particle*> all_truth_particles;
  std::map<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster>> alltruthclusters;
  for(auto trkid:trkrlist){

  for (const auto &hitsetkey :
       clustermap->getHitSetKeys(trkid.second)) {
    auto range = clustermap->getClusters(hitsetkey);
    for (auto clusterIter = range.first; clusterIter != range.second; ++clusterIter){
      const auto &key = clusterIter->first;
      const auto &cluster = clusterIter->second;

      //uint8_t layer = TrkrDefs::getLayer(key);
      //std::cout<<"layer: "<<(int)layer<<std::endl;

      std::set<PHG4Particle*> truth_withcluster = clustereval->all_truth_particles(key);
      //check if the set it empty
      if(truth_withcluster.empty()){
        std::cout << "calotrkana::process_event(PHCompositeNode *topNode) truth_withcluster is empty" << std::endl;
      }
      //merge the two sets
      all_truth_particles.insert(truth_withcluster.begin(), truth_withcluster.end());


      const auto global = m_tGeometry->getGlobalPosition(key, cluster);

      float x = global(0);
      float y = global(1);
      float z = global(2);

      float e = cluster->getAdc();
      //float t = cluster->getTime();
      float t = 0;
      //strobe for mvtx, crossing for intt
      if(trkid.second == TrkrDefs::TrkrId::mvtxId){
        t = MvtxDefs::getStrobeId(key);
      }
      else if(trkid.second == TrkrDefs::TrkrId::inttId){
        t = InttDefs::getTimeBucketId(key); 
      }
      else{
        t = 0;
      }
      // find the truth cluster
      std::pair<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster>> truth_cluster_pair = clustereval->max_truth_cluster_by_energy(key);
      TrkrDefs::cluskey truth_ckey = truth_cluster_pair.first;
      std::shared_ptr<TrkrCluster> truth_cluster = truth_cluster_pair.second;
      //put it into the truth cluster map shift the id up by 1 and left 0 for unmatched clusters
      unsigned int clusterid = (unsigned int)TrkrDefs::getClusIndex(truth_ckey);
      //check if truth cluster is nullptr
      if(truth_cluster){
        alltruthclusters[truth_ckey] = truth_cluster;
        clusterid++;
      }
      /*
      else{
        std::cout << "calotrkana::process_event(PHCompositeNode *topNode) truth cluster is nullptr" << std::endl;
      }
      */
      m_reco_cluster_E[m_nRecoClusters] = e;
      m_reco_cluster_x[m_nRecoClusters] = x;
      m_reco_cluster_y[m_nRecoClusters] = y;
      m_reco_cluster_z[m_nRecoClusters] = z;
      m_reco_cluster_t[m_nRecoClusters] = t;
      m_reco_cluster_detid[m_nRecoClusters] = static_cast<int>(TrkrDefs::getTrkrId(key));
      m_reco_cluster_trcluster_id[m_nRecoClusters] = clusterid;
      m_nRecoClusters++;

      if(m_nRecoClusters >= trackrecoclustermaxlength){
        std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nRecoClusters exceeds max length" << std::endl;
        exit(1);
      }
     
    }
  }
}
//loop over all truth particles and find all associated truth clusters(just in case there are truth cluster that are not associated by any reco clusters)
//maybe unnecessary...just to be safe :)
for(auto truth:all_truth_particles){
  std::map<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster>> truth_clusters = trutheval->all_truth_clusters(truth);
  //merge the two maps
  alltruthclusters.insert(truth_clusters.begin(), truth_clusters.end());
}

for(auto truth:alltruthclusters){
  const auto &key = truth.first;
  const auto &cluster = truth.second;

  float x = cluster->getPosition(0);
  float y = cluster->getPosition(1);
  float z = cluster->getPosition(2);

  float e = cluster->getAdc();
  //float t = cluster->getTime();
  float t = 0;
  //strobe for mvtx, crossing for intt,(however this is all zero for truth clusters?)
  //https://github.com/sPHENIX-Collaboration/coresoftware/blob/master/simulation/g4simulation/g4eval/SvtxTruthEval.cc#L351
  if(TrkrDefs::getTrkrId(key) == TrkrDefs::TrkrId::mvtxId){
    t = MvtxDefs::getStrobeId(key);
  }
  else if(TrkrDefs::getTrkrId(key) == TrkrDefs::TrkrId::inttId){
    t = InttDefs::getTimeBucketId(key); 
  }
  else{
    t = 0;
  }
  //get truth G4Hit
  std::set<PHG4Hit*> g4hits = trutheval->get_truth_hits_from_truth_cluster(key);
  //for truth clusters all hits are belong to the same particle
  //check if g4hits is empty
  if(g4hits.empty()){
    std::cout << "calotrkana::process_event(PHCompositeNode *topNode) g4hits for a truth cluster is empty something is very wrong" << std::endl;
    exit(1);
  }

  PHG4Particle *truth_particle = trutheval->get_particle(*g4hits.begin());
  int truth_particle_trackid = truth_particle->get_track_id();
  //std::cout<<"truth_particle_trackid: "<<truth_particle_trackid<<std::endl;
  
  //cluster shift up by 1
  unsigned int clusterid = (unsigned int)TrkrDefs::getClusIndex(key) + 1;
  m_truth_cluster_E[m_nTruthClusters] = e;
  m_truth_cluster_x[m_nTruthClusters] = x;
  m_truth_cluster_y[m_nTruthClusters] = y;
  m_truth_cluster_z[m_nTruthClusters] = z;
  m_truth_cluster_t[m_nTruthClusters] = t;
  m_truth_cluster_id[m_nTruthClusters] = clusterid;
  m_truth_cluster_trparticle_track_id[m_nTruthClusters] = truth_particle_trackid;
  //std::cout<<"truth_particle_trackid: "<<truth_particle_trackid<<std::endl;
  m_truth_cluster_detid[m_nTruthClusters] = static_cast<int>(TrkrDefs::getTrkrId(key));
  m_nTruthClusters++;

  if(m_nTruthClusters >= truthclustermaxlength){
    std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nTruthClusters exceeds max length" << std::endl;
    exit(1);
  }
}




//get (secondary) truth particles associated with clusters
all_truth_particles.insert(primary_particles.begin(), primary_particles.end());
for(auto truth:all_truth_particles){
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
  m_particle_track_id[m_nParticles] = truth->get_track_id();
  m_particle_primary_id[m_nParticles] = truth->get_primary_id();
  m_nParticles++;

  if(m_nParticles >= ptruthmaxlength){
    std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nParticles exceeds max length" << std::endl;
    exit(1);
  }
}

  T->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int calotrkana::ResetEvent(PHCompositeNode *topNode) {

  m_nHits = 0;
  m_nParticles = 0;
  m_nRecoClusters = 0;
  m_nTruthClusters = 0;

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
