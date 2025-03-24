//____________________________________________________________________________..
//
//
//____________________________________________________________________________..

#include "calotrkana.h"

// Fun4all stuff
#include <ffaobjects/EventHeader.h>
#include <ffaobjects/RunHeader.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

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
#include <TObject.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

// just in case in the future I want to get HEPMC info
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

#include <centrality/CentralityInfo.h>

#include <epd/EpdGeom.h>

// tracking stuff
#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h> // for getTrkrId, getHit...
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/InttDefs.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <g4eval/SvtxEvalStack.h>

//____________________________________________________________________________..
calotrkana::calotrkana(const std::string &name = "calotrkana",
                       const std::string &outName = "ZDCresult")
    : SubsysReco(name), Outfile(outName)
{
  std::cout << "calotrkana::calotrkana(const std::string &name) Calling ctor"
            << std::endl;
}

//____________________________________________________________________________..
calotrkana::~calotrkana()
{
  std::cout << "calotrkana::~calotrkana() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int calotrkana::Init(PHCompositeNode *topNode)
{

  std::cout << "calotrkana::Init(PHCompositeNode *topNode) Initializing"
            << std::endl;
  out = new TFile(Outfile.c_str(), "RECREATE");
  T = new TTree("T", "Tree for Reco Hits and truth");

  // global branches
  // HEPMC HI info
  T->Branch("b", &m_b, "b/F");
  T->Branch("b_phi", &m_b_phi, "b_phi/F");
  T->Branch("Ncoll", &m_Ncoll, "Ncoll/I");
  T->Branch("Ncoll_hard", &m_Ncoll_hard, "Ncoll_hard/I");
  T->Branch("Npart_proj", &m_Npart_proj, "Npart_proj/I");
  T->Branch("Npart_targ", &m_Npart_targ, "Npart_targ/I");
  // centrality
  T->Branch("centile", &m_cent, "centile/F");
  // event meta data
  T->Branch("runnumber", &m_runnumber, "runnumber/I");
  T->Branch("evtnumber", &m_evtnumber, "evtnumber/I");
  // string for filename
  T->Branch("filename", &m_filename);

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
  T->Branch("particle_parent_id", &m_particle_parent_id, "particle_parent_id[nParticles]/I");
  T->Branch("particle_is_pythia_primary", &m_particle_is_pythia_primary, "particle_is_pythia_primary[nParticles]/I");
  T->Branch("particle_is_embedded", &m_particle_is_embedded, "particle_is_embedded[nParticles]/I");
  T->Branch("particle_charge", &m_particle_charge, "particle_charge[nParticles]/I");
  T->Branch("nRecoClusters", &m_nRecoClusters, "nRecoClusters/I");
  T->Branch("reco_cluster_E", &m_reco_cluster_E, "reco_cluster_E[nRecoClusters]/F");
  T->Branch("reco_cluster_x", &m_reco_cluster_x, "reco_cluster_x[nRecoClusters]/F");
  T->Branch("reco_cluster_y", &m_reco_cluster_y, "reco_cluster_y[nRecoClusters]/F");
  T->Branch("reco_cluster_z", &m_reco_cluster_z, "reco_cluster_z[nRecoClusters]/F");
  T->Branch("reco_cluster_t", &m_reco_cluster_t, "reco_cluster_t[nRecoClusters]/F");
  T->Branch("reco_cluster_detid", &m_reco_cluster_detid, "reco_cluster_detid[nRecoClusters]/I");
  T->Branch("reco_cluster_id", &m_reco_cluster_id, "reco_cluster_id[nRecoClusters]/l");
  T->Branch("reco_cluster_trcluster_id", &m_reco_cluster_trcluster_id, "reco_cluster_trcluster_id[nRecoClusters]/i");
  T->Branch("reco_cluster_g4hit_id", &m_reco_cluster_g4hit_id, "reco_cluster_g4hit_id[nRecoClusters]/l");
  T->Branch("nTruthClusters", &m_nTruthClusters, "nTruthClusters/I");
  T->Branch("truth_cluster_E", &m_truth_cluster_E, "truth_cluster_E[nTruthClusters]/F");
  T->Branch("truth_cluster_x", &m_truth_cluster_x, "truth_cluster_x[nTruthClusters]/F");
  T->Branch("truth_cluster_y", &m_truth_cluster_y, "truth_cluster_y[nTruthClusters]/F");
  T->Branch("truth_cluster_z", &m_truth_cluster_z, "truth_cluster_z[nTruthClusters]/F");
  T->Branch("truth_cluster_t", &m_truth_cluster_t, "truth_cluster_t[nTruthClusters]/F");
  T->Branch("truth_cluster_detid", &m_truth_cluster_detid, "truth_cluster_detid[nTruthClusters]/I");
  T->Branch("truth_cluster_id", &m_truth_cluster_id, "truth_cluster_id[nTruthClusters]/i");
  T->Branch("truth_cluster_trparticle_track_id", &m_truth_cluster_trparticle_track_id, "truth_cluster_trparticle_track_id[nTruthClusters]/I");
  T->Branch("nTrackG4Hits", &m_nTrackG4Hits, "nTrackG4Hits/I");
  T->Branch("track_g4hit_x", &m_track_g4hit_x, "track_g4hit_x[nTrackG4Hits]/F");
  T->Branch("track_g4hit_y", &m_track_g4hit_y, "track_g4hit_y[nTrackG4Hits]/F");
  T->Branch("track_g4hit_z", &m_track_g4hit_z, "track_g4hit_z[nTrackG4Hits]/F");
  T->Branch("track_g4hit_t", &m_track_g4hit_t, "track_g4hit_t[nTrackG4Hits]/F");
  T->Branch("track_g4hit_E", &m_track_g4hit_E, "track_g4hit_E[nTrackG4Hits]/F");
  T->Branch("track_g4hit_trparticle_track_id", &m_track_g4hit_trparticle_track_id, "track_g4hit_trparticle_track_id[nTrackG4Hits]/I");
  T->Branch("track_g4hit_id", &m_track_g4hit_id, "track_g4hit_id[nTrackG4Hits]/l");

  // TPC seeds
  /*
  T->Branch("tpc_seeds_id", &m_tpc_seeds_id);
  T->Branch("tpc_seeds_nclusters", &m_tpc_seeds_nclusters);
  T->Branch("tpc_seeds_start_idx", &m_tpc_seeds_start_idx);
  T->Branch("tpc_seeds_clusters", &m_tpc_seeds_clusters);
  */
  T->Branch("nTPCSeeds", &m_nTPCSeeds, "nTPCSeeds/I");
  T->Branch("nTPCSeedsClusters", &m_nTPCSeedsClusters, "nTPCSeedsClusters/I");
  T->Branch("tpc_seeds_id", &m_tpc_seeds_id, "tpc_seeds_id[nTPCSeeds]/i");
  T->Branch("tpc_seeds_nclusters", &m_tpc_seeds_nclusters, "tpc_seeds_nclusters[nTPCSeeds]/i");
  T->Branch("tpc_seeds_start_idx", &m_tpc_seeds_start_idx, "tpc_seeds_start_idx[nTPCSeeds]/i");
  T->Branch("tpc_seeds_clusters", &m_tpc_seeds_clusters, "tpc_seeds_clusters[nTPCSeedsClusters]/l");

  _pdg = new TDatabasePDG();
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Print("NODETREE");

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int calotrkana::InitRun(PHCompositeNode *topNode)
{
  std::cout
      << "calotrkana::InitRun(PHCompositeNode *topNode) Initializing for Run XXX"
      << std::endl;

  // get runnumber here
  RunHeader *runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (!runheader)
  {
    std::cout << "can't find runheader" << std::endl;
    return 1;
  }
  m_runnumber = runheader->get_RunNumber();
  // get the filename flag via rc
  //!!!!!!!!!!!!!! this only works if we read in one segment at a time...I don't think there is a simple way to get the segment name in F4A...
  recoConsts *rc = recoConsts::instance();
  m_filename = rc->get_StringFlag("Sim_File_Name");
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int calotrkana::process_event(PHCompositeNode *topNode)
{

  if (m_HI)
  {
    // hepmc info
    PHHepMCGenEventMap *genevtmap =
        findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
    if (!genevtmap)
    {
      std::cout << "no genevtmap" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    for (PHHepMCGenEventMap::Iter iter = genevtmap->begin();
         iter != genevtmap->end(); ++iter)
    {
      PHHepMCGenEvent *genevt = iter->second;
      // check if embedded
      if (genevt->get_embedding_id() != 0)
        continue;
      HepMC::GenEvent *event = genevt->getEvent();
      if (!event)
      {
        std::cout << PHWHERE << " no evt pointer under HEPMC Node found"
                  << std::endl;
      }
      else
      {
        HepMC::HeavyIon *hi = event->heavy_ion();
        if (!hi)
        {
          std::cout << PHWHERE
                    << ": no heavy ion info found in hepmc event, skip this event"
                    << std::endl;
        }
        m_b = hi->impact_parameter();
        m_b_phi = hi->event_plane_angle();
        m_Ncoll = hi->Ncoll();
        m_Ncoll_hard = hi->Ncoll_hard();
        m_Npart_proj = hi->Npart_proj();
        m_Npart_targ = hi->Npart_targ();
        std::cout<<"npart_proj: "<<m_Npart_proj<<" npart_targ: "<<m_Npart_targ<<std::endl;
        int Npart_total = m_Npart_proj + m_Npart_targ;
        if(Npart_total> m_Npart_max || Npart_total < m_Npart_min)
        {
          std::cout<<"Npart_total: "<<Npart_total<<" is out of range: "<<m_Npart_min<<" - "<<m_Npart_max<<std::endl;
          std::cout<<"Skipping event"<<std::endl;
          return Fun4AllReturnCodes::EVENT_OK;
        }
      }
    }
    // centrality and vertex maybe
    CentralityInfo *centrality_info = findNode::getClass<CentralityInfo>(topNode, "Centrality");
    if (!centrality_info)
    {
      std::cout << "calotrkana::process_event(PHCompositeNode *topNode) No centrality info found"
                << std::endl;
    }
    else
    {
      // get cent
      m_cent = (centrality_info->has_centile(CentralityInfo::PROP::mbd_NS) ? centrality_info->get_centile(CentralityInfo::PROP::mbd_NS) : -999.99);
    }
  }

  // truth container here
  PHG4TruthInfoContainer *truthinfo =
      findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truthinfo)
  {
    std::cout << "calotrkana::process_event(PHCompositeNode *topNode) No truth "
                 "info found"
              << std::endl;
  }

  // loop over truth primary particles
  PHG4TruthInfoContainer::ConstRange range =
      truthinfo->GetPrimaryParticleRange();
  // set of primary particles
  std::set<PHG4Particle *> primary_particles;
  for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first;
       truth_itr != range.second; ++truth_itr)
  {
    PHG4Particle *truth = truth_itr->second;
    primary_particles.insert(truth);
  }
  // loop over secondary particles
  PHG4TruthInfoContainer::ConstRange range_sec =
      truthinfo->GetSecondaryParticleRange();
  // set of secondary particles
  std::set<PHG4Particle *> secondary_particles;
  float max_r = 93;
  float max_z = 150;
  for (PHG4TruthInfoContainer::ConstIterator truth_itr = range_sec.first;
       truth_itr != range_sec.second; ++truth_itr)
  {
    PHG4Particle *truth = truth_itr->second;
    int vtxid = truth->get_vtx_id();
    PHG4VtxPoint *vtx = truthinfo->GetVtx(vtxid);
    float vtx_x = vtx->get_x();
    float vtx_y = vtx->get_y();
    float vtx_z = vtx->get_z();
    float vtx_r = sqrt(vtx_x * vtx_x + vtx_y * vtx_y);
    if (vtx_r > max_r)
    {
      continue;
    }
    if (abs(vtx_z) > max_z)
    {
      continue;
    }
    secondary_particles.insert(truth);
  }

  // CEMC
  TowerInfoContainer *CEMC_towers_sim =
      findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  RawTowerGeomContainer *CEMC_geom =
      findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");

  int CEMCsize = CEMC_towers_sim->size();
  for (int i = 0; i < CEMCsize; i++)
  {
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
    m_Hit_t[m_nHits] = tower->get_time_float() * sampletons;
    m_Hit_detid[m_nHits] = static_cast<int>(calotrkana::cemcId);
    m_nHits++;

    if (m_nHits >= recomaxlength)
    {
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
  for (int i = 0; i < IHCALsize; i++)
  {
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
    m_Hit_t[m_nHits] = tower->get_time_float() * sampletons;
    m_Hit_detid[m_nHits] = static_cast<int>(calotrkana::ihcalId);
    m_nHits++;

    if (m_nHits >= recomaxlength)
    {
      std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nHits exceeds max length" << std::endl;
      exit(1);
    }
  }

  // OHCAL
  TowerInfoContainer *OHCAL_towers_sim =
      findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  RawTowerGeomContainer *OHCAL_geom =
      findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  if (!OHCAL_towers_sim || !OHCAL_geom)
  {
    std::cout << "calotrkana::process_event(PHCompositeNode *topNode) No OHCAL "
                 "towers found"
              << std::endl;
  }
  int OHCALsize = OHCAL_towers_sim->size();
  for (int i = 0; i < OHCALsize; i++)
  {
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
    m_Hit_t[m_nHits] = tower->get_time_float() * sampletons;
    m_Hit_detid[m_nHits] = static_cast<int>(calotrkana::ohcalId);
    m_nHits++;

    if (m_nHits >= recomaxlength)
    {
      std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nHits exceeds max length" << std::endl;
      exit(1);
    }
  }

  // SEPD
  TowerInfoContainer *EPD_towers_sim =
      findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_EPD");
  EpdGeom *EPD_geom =
      findNode::getClass<EpdGeom>(topNode, "TOWERGEOM_EPD");
  if (!OHCAL_towers_sim || !OHCAL_geom)
  {
    std::cout << "calotrkana::process_event(PHCompositeNode *topNode) No EPD "
                 "towers found"
              << std::endl;
  }

  unsigned int EPDsize = EPD_towers_sim->size();
  for (unsigned int i = 0; i < EPDsize; i++)
  {
    TowerInfo *tower = EPD_towers_sim->get_tower_at_channel(i);
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

    if (m_nHits >= recomaxlength)
    {
      std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nHits exceeds max length" << std::endl;
      exit(1);
    }
  }

  // MBD
  MbdPmtContainer *mbdpmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  MbdGeom *mbdgeom = findNode::getClass<MbdGeom>(topNode, "MbdGeom");
  if (!mbdpmts || !mbdgeom)
  {
    std::cout << "calotrkana::process_event(PHCompositeNode *topNode) No MBD "
                 "pmts found"
              << std::endl;
  }
  for (int ipmt = 0; ipmt < mbdpmts->get_npmt(); ipmt++)
  {
    float mbdq = mbdpmts->get_pmt(ipmt)->get_q();
    if (mbdq < 0.25) // from MBD reco hit threshold
      continue;
    float mbdt = mbdpmts->get_pmt(ipmt)->get_time();
    // check if mbdt is undefined
    if (mbdt > 1E6)
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

    if (m_nHits >= recomaxlength)
    {
      std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nHits exceeds max length" << std::endl;
      exit(1);
    }
  }
  // find secondary particles associated with tack clusters
  if (!m_svtxEvalStack)
  {
    m_svtxEvalStack = new SvtxEvalStack(topNode);
  }
  // do eval
  m_svtxEvalStack->next_event(topNode);

  /// Get the cluster evaluator
  SvtxClusterEval *clustereval = m_svtxEvalStack->get_cluster_eval();
  // clustereval->set_verbosity(1);
  SvtxTruthEval *trutheval = m_svtxEvalStack->get_truth_eval();
  // trutheval->set_verbosity(2);
  //  tracking
  ActsGeometry *m_tGeometry =
      findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE << "No acts tracking geometry, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  TrkrClusterContainer *clustermap =
      findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!clustermap)
  {
    std::cout << "calotrkana::process_event(PHCompositeNode *topNode) No cluster "
                 "map found"
              << std::endl;
  }
  /*
  TrkrClusterContainer *truthclustermap =
      findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_TRUTHCLUSTERCONTAINER");
  if (!truthclustermap) {
    std::cout << "calotrkana::process_event(PHCompositeNode *topNode) No truth cluster "
                 "map found"
              << std::endl;
  }
  */
  std::set<PHG4Particle *> all_truth_particles;
  std::set<PHG4Hit *> all_track_g4hits;
  std::map<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster>> alltruthclusters;
  std::cout << "calotrkana::process_event(PHCompositeNode *topNode) clustermap size: " << clustermap->size() << std::endl;
  for (auto trkid : trkrlist)
  {

    for (const auto &hitsetkey :
         clustermap->getHitSetKeys(trkid.second))
    {
      auto range = clustermap->getClusters(hitsetkey);
      // this part somehow currpt the memory...
      /*
      auto truthrange = truthclustermap->getClusters(hitsetkey);
      for (auto clusterIter = truthrange.first; clusterIter != truthrange.second; ++clusterIter){
        const auto &key = clusterIter->first;
        //const auto &cluster = clusterIter->second;
        std::shared_ptr<TrkrCluster> truth_cluster = std::make_shared<TrkrCluster>(*clusterIter->second);
        //make pair
        alltruthclusters.insert(std::make_pair(key, truth_cluster));
      }
      */
      for (auto clusterIter = range.first; clusterIter != range.second; ++clusterIter)
      {
        const auto &key = clusterIter->first;
        const auto &cluster = clusterIter->second;

        // uint8_t layer = TrkrDefs::getLayer(key);
        // std::cout<<"layer: "<<(int)layer<<std::endl;

        std::set<PHG4Particle *> truth_withcluster = clustereval->all_truth_particles(key);
        // check if the set it empty
        if (truth_withcluster.empty())
        {
          std::cout << "calotrkana::process_event(PHCompositeNode *topNode) truth_withcluster is empty" << std::endl;
        }
        // merge the two sets
        all_truth_particles.insert(truth_withcluster.begin(), truth_withcluster.end());

        std::set<PHG4Hit *> g4hits = clustereval->all_truth_hits(key);
        if (g4hits.empty())
        {
          std::cout << "calotrkana::process_event(PHCompositeNode *topNode) g4hits is empty" << std::endl;
        }
        all_track_g4hits.insert(g4hits.begin(), g4hits.end());

        const auto global = m_tGeometry->getGlobalPosition(key, cluster);

        float x = global(0);
        float y = global(1);
        float z = global(2);
        /*
        int truthtrackid = (*truth_withcluster.begin())->get_track_id();
        if(abs(z)>200 || truthtrackid == -1134){
          std::cout<<"x: "<<x<<" y: "<<y<<" z: "<<z<<std::endl;

          std::cout<<"truthtrackid: "<<truthtrackid<<std::endl;
          float localy = cluster->getLocalY();
          std::cout<<"localy: "<<localy<<std::endl;
        }
        */

        float e = cluster->getAdc();
        // float t = cluster->getTime();
        float t = 0;
        // strobe for mvtx, crossing for intt
        if (trkid.second == TrkrDefs::TrkrId::mvtxId)
        {
          t = MvtxDefs::getStrobeId(key);
        }
        else if (trkid.second == TrkrDefs::TrkrId::inttId)
        {
          t = InttDefs::getTimeBucketId(key);
        }
        else
        {
          t = 0;
        }
        // find the truth cluster
        std::pair<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster>> truth_cluster_pair = clustereval->max_truth_cluster_by_energy(key);
        TrkrDefs::cluskey truth_ckey = truth_cluster_pair.first;
        std::shared_ptr<TrkrCluster> truth_cluster = truth_cluster_pair.second;
        // put it into the truth cluster map shift the id up by 1 and left 0 for unmatched clusters
        unsigned int clusterid = (unsigned int)TrkrDefs::getClusIndex(truth_ckey);
        // check if truth cluster is nullptr
        if (truth_cluster)
        {
          alltruthclusters[truth_ckey] = truth_cluster;
          clusterid++;
        }
        else
        {
          // not matched set to NULL
          clusterid = 0;
        }
        // find the best matched G4Hit
        PHG4Hit *g4hit = clustereval->max_truth_hit_by_energy(key);
        ULong64_t g4hit_id = 0;
        if (g4hit)
        {
          g4hit_id = g4hit->get_hit_id();
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
        // key here is uint64_t (https://github.com/sPHENIX-Collaboration/coresoftware/blob/master/offline/packages/trackbase/TrkrDefs.h#L27)
        m_reco_cluster_id[m_nRecoClusters] = static_cast<ULong64_t>(key);
        m_reco_cluster_trcluster_id[m_nRecoClusters] = clusterid;
        m_reco_cluster_g4hit_id[m_nRecoClusters] = g4hit_id;
        m_nRecoClusters++;

        if (m_nRecoClusters >= trackrecoclustermaxlength)
        {
          std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nRecoClusters exceeds max length" << std::endl;
          exit(1);
        }
      }
    }
  }

  // loop over all TPC seeds
  /*
  m_tpc_seeds_id.clear();
  m_tpc_seeds_nclusters.clear();
  m_tpc_seeds_start_idx.clear();
  m_tpc_seeds_clusters.clear();
  int tpc_seed_start_idx = 0;
  SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  for (const auto &[key, track] : *trackmap)
  {
    if (!track)
    {
      continue;
    }
    std::vector<ULong64_t> cluster_keys;
    int seed_ntpc_clusters = 0;

    for (const auto &ckey : get_cluster_keys(track))
    {
      if (TrkrDefs::getTrkrId(ckey) == TrkrDefs::tpcId)
      {
        cluster_keys.push_back(static_cast<ULong64_t>(ckey));
        m_tpc_seeds_clusters.push_back(static_cast<ULong64_t>(ckey));
        seed_ntpc_clusters++;
      }
    }
    if (seed_ntpc_clusters == 0)
    {
      continue;
    }

    m_tpc_seeds_id.push_back(static_cast<unsigned int>(key));
    m_tpc_seeds_nclusters.push_back(seed_ntpc_clusters);
    m_tpc_seeds_start_idx.push_back(tpc_seed_start_idx);
    tpc_seed_start_idx += seed_ntpc_clusters;
    //m_tpc_seeds_clusters.push_back(cluster_keys);
  }
  */

  // loop over all TPC seeds (use raw arrrays)
  int tpc_seed_start_idx = 0;
  SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  for (const auto &[key, track] : *trackmap)
  {
    if (!track)
    {
      continue;
    }
    int seed_ntpc_clusters = 0;

    for (const auto &ckey : get_cluster_keys(track))
    {
      if (TrkrDefs::getTrkrId(ckey) == TrkrDefs::tpcId)
      {

        m_tpc_seeds_clusters[m_nTPCSeedsClusters] = static_cast<ULong64_t>(ckey);
        m_nTPCSeedsClusters++;
        seed_ntpc_clusters++;
        if (m_nTPCSeedsClusters >= tpcseedclustermaxlength)
        {
          std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nTPCSeedsClusters exceeds max length" << std::endl;
          exit(1);
        }
      }
    }
    if (seed_ntpc_clusters == 0)
    {
      continue;
    }

    m_tpc_seeds_id[m_nTPCSeeds] = static_cast<unsigned int>(key);
    m_tpc_seeds_nclusters[m_nTPCSeeds] = seed_ntpc_clusters;
    m_tpc_seeds_start_idx[m_nTPCSeeds] = tpc_seed_start_idx;
    m_nTPCSeeds++;
    tpc_seed_start_idx += seed_ntpc_clusters;

    if (m_nTPCSeeds >= tpcseedmaxlength)
    {
      std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nTPCSeeds exceeds max length" << std::endl;
      exit(1);
    }
  }
  // loop over all truth particles and find all associated truth clusters(just in case there are truth cluster that are not associated by any reco clusters)
  // maybe unnecessary...just to be safe :)
  for (auto truth : all_truth_particles)
  {
    std::map<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster>> truth_clusters = trutheval->all_truth_clusters(truth);
    // merge the two maps
    alltruthclusters.insert(truth_clusters.begin(), truth_clusters.end());
  }

  for (auto truth : alltruthclusters)
  {
    const auto &key = truth.first;
    const auto &cluster = truth.second;
    if (!cluster)
    {
      std::cout << "calotrkana::process_event(PHCompositeNode *topNode) cluster is nullptr" << std::endl;
      continue;
    }
    float x = cluster->getPosition(0);
    float y = cluster->getPosition(1);
    float z = cluster->getPosition(2);

    float e = cluster->getAdc();
    // float t = cluster->getTime();
    float t = 0;
    // strobe for mvtx, crossing for intt,(however this is all zero for truth clusters?)
    // https://github.com/sPHENIX-Collaboration/coresoftware/blob/master/simulation/g4simulation/g4eval/SvtxTruthEval.cc#L351
    if (TrkrDefs::getTrkrId(key) == TrkrDefs::TrkrId::mvtxId)
    {
      t = MvtxDefs::getStrobeId(key);
    }
    else if (TrkrDefs::getTrkrId(key) == TrkrDefs::TrkrId::inttId)
    {
      t = InttDefs::getTimeBucketId(key);
    }
    else
    {
      t = 0;
    }
    // get truth G4Hit
    std::set<PHG4Hit *> g4hits = trutheval->get_truth_hits_from_truth_cluster(key);
    // for truth clusters all hits are belong to the same particle
    // check if g4hits is empty
    int truth_particle_trackid = 0;
    if (!g4hits.empty())
    {
      PHG4Particle *truth_particle = trutheval->get_particle(*g4hits.begin());
      truth_particle_trackid = truth_particle->get_track_id();
    }
    else
    {
      continue;
    }

    // std::cout<<"truth_particle_trackid: "<<truth_particle_trackid<<std::endl;

    // cluster shift up by 1
    unsigned int clusterid = (unsigned int)TrkrDefs::getClusIndex(key) + 1;
    m_truth_cluster_E[m_nTruthClusters] = e;
    m_truth_cluster_x[m_nTruthClusters] = x;
    m_truth_cluster_y[m_nTruthClusters] = y;
    m_truth_cluster_z[m_nTruthClusters] = z;
    m_truth_cluster_t[m_nTruthClusters] = t;
    m_truth_cluster_id[m_nTruthClusters] = clusterid;
    m_truth_cluster_trparticle_track_id[m_nTruthClusters] = truth_particle_trackid;
    // std::cout<<"truth_particle_trackid: "<<truth_particle_trackid<<std::endl;
    m_truth_cluster_detid[m_nTruthClusters] = static_cast<int>(TrkrDefs::getTrkrId(key));
    m_nTruthClusters++;

    if (m_nTruthClusters >= truthclustermaxlength)
    {
      std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nTruthClusters exceeds max length" << std::endl;
      exit(1);
    }
  }
  std::set<PHG4Hit *> alltruthhits = trutheval->all_truth_hits();

  all_track_g4hits.insert(alltruthhits.begin(), alltruthhits.end());

  // loop over all G4Hits associated with track clusters
  for (auto g4hit : all_track_g4hits)
  {
    float x = g4hit->get_avg_x();
    float y = g4hit->get_avg_y();
    float z = g4hit->get_avg_z();
    float e = g4hit->get_edep();
    float t = g4hit->get_avg_t();
    PHG4Particle *truth_particle = trutheval->get_particle(g4hit);
    all_truth_particles.insert(truth_particle);
    int trparticle_track_id = truth_particle->get_track_id();
    m_track_g4hit_x[m_nTrackG4Hits] = x;
    m_track_g4hit_y[m_nTrackG4Hits] = y;
    m_track_g4hit_z[m_nTrackG4Hits] = z;
    m_track_g4hit_t[m_nTrackG4Hits] = t;
    m_track_g4hit_E[m_nTrackG4Hits] = e;
    m_track_g4hit_trparticle_track_id[m_nTrackG4Hits] = trparticle_track_id;
    m_track_g4hit_id[m_nTrackG4Hits] = g4hit->get_hit_id();
    m_nTrackG4Hits++;

    if (m_nTrackG4Hits >= truthtrackg4hitmaxlength)
    {
      std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nTrackG4Hits exceeds max length" << std::endl;
      exit(1);
    }
  }
  // extra step to get all secondaries parent history upstream to the saved particles
  std::set<PHG4Particle *> completed_particle_tree;
  all_truth_particles.insert(secondary_particles.begin(), secondary_particles.end());
  for (auto truth : all_truth_particles)
  {
    int parent_id = truth->get_parent_id();
    PHG4Particle *parent = truthinfo->GetParticle(parent_id);
    while (parent)
    {
      // if parent is already in the set, no need to go further
      if (completed_particle_tree.find(parent) != completed_particle_tree.end())
      {
        break;
      }
      completed_particle_tree.insert(parent);
      parent_id = parent->get_parent_id();
      parent = truthinfo->GetParticle(parent_id);
    }
  }

  // get (secondary) truth particles associated with clusters
  all_truth_particles.insert(completed_particle_tree.begin(), completed_particle_tree.end());
  all_truth_particles.insert(primary_particles.begin(), primary_particles.end());

  for (auto truth : all_truth_particles)
  {

    int is_pythia_primary = trutheval->is_primary(truth);
    int is_embedded = trutheval->get_embed(truth);
    int pid = truth->get_pid();
    TParticlePDG *partinfo = _pdg->GetParticle(pid);
    float charge = std::nan("");
    if (partinfo)
    {
      charge = partinfo->Charge() / 3; // PDG gives charge in 1/3 e
    }
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
    m_particle_parent_id[m_nParticles] = truth->get_parent_id();
    m_particle_is_pythia_primary[m_nParticles] = is_pythia_primary;
    m_particle_is_embedded[m_nParticles] = is_embedded;
    m_particle_charge[m_nParticles] = (int)charge;
    m_nParticles++;

    if (m_nParticles >= ptruthmaxlength)
    {
      std::cout << "calotrkana::process_event(PHCompositeNode *topNode) m_nParticles exceeds max length" << std::endl;
      exit(1);
    }
  }

  T->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int calotrkana::ResetEvent(PHCompositeNode *topNode)
{

  m_nHits = 0;
  m_nParticles = 0;
  m_nRecoClusters = 0;
  m_nTruthClusters = 0;
  m_nTrackG4Hits = 0;
  m_nTPCSeeds = 0;
  m_nTPCSeedsClusters = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int calotrkana::End(PHCompositeNode *topNode)
{
  out->cd();

  // T->Write();
  T->Write("", TObject::kOverwrite);

  out->Close();
  delete out;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int calotrkana::Reset(PHCompositeNode *topNode)
{
  std::cout << "calotrkana::Reset(PHCompositeNode *topNode) being Reset"
            << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void calotrkana::Print(const std::string &what) const
{
  std::cout << "calotrkana::Print(const std::string &what) const Printing info for "
            << what << std::endl;
}

std::vector<TrkrDefs::cluskey> calotrkana::get_cluster_keys(SvtxTrack *track)
{
  std::vector<TrkrDefs::cluskey> out;
  for (const auto &seed : {track->get_silicon_seed(), track->get_tpc_seed()})
  {
    if (seed)
    {
      std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
    }
  }
  return out;
}