// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef calotrkana_H
#define calotrkana_H

#include <fun4all/SubsysReco.h>
// ROOT stuff
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <calobase/TowerInfoDefs.h>
#include <caloreco/CaloTowerDefs.h>

#include <trackbase/TrkrDefs.h>

#include <string>

namespace HepMC {
class GenEvent;
}

class PHCompositeNode;
class Fun4AllHistoManager;
class PHG4Hit;
class PHG4CylinderCellGeom_Spacalv1;
class PHG4CylinderGeom_Spacalv3;
class SvtxEvalStack;
class TDatabasePDG;

class calotrkana : public SubsysReco {
public:
  calotrkana(const std::string &name, const std::string &outName);

  ~calotrkana() override;

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;


private:
  TTree *T = nullptr;
  const float sampletons = 50./3.;
  // all reco stuff
  static const int recomaxlength = 1E7;
  
  float m_Hit_E[recomaxlength] = {0};
  float m_Hit_x[recomaxlength] = {0};
  float m_Hit_y[recomaxlength] = {0};
  float m_Hit_z[recomaxlength] = {0};
  float m_Hit_t[recomaxlength] = {0};
  int m_Hit_detid[recomaxlength] = {0};
  int m_nHits = 0;
  //I'm seperating the traking clusters from the rest
  static const int trackrecoclustermaxlength = 1E7;
  float m_reco_cluster_E[trackrecoclustermaxlength] = {0};
  float m_reco_cluster_x[trackrecoclustermaxlength] = {0};
  float m_reco_cluster_y[trackrecoclustermaxlength] = {0};
  float m_reco_cluster_z[trackrecoclustermaxlength] = {0};
  //only INTT has a crossing id and MVTX has a strobe id(but for no pileup it's probably always 0)
  float m_reco_cluster_t[trackrecoclustermaxlength] = {0};
  int m_reco_cluster_detid[trackrecoclustermaxlength] = {0};
  //best matching truth cluster
  unsigned int m_reco_cluster_trcluster_id[trackrecoclustermaxlength] = {0};
  //best matching G4Hit
  ULong64_t m_reco_cluster_g4hit_id[trackrecoclustermaxlength] = {0};
  int m_nRecoClusters = 0;
  //truth cluster
  static const int truthclustermaxlength = 1E7;
  float m_truth_cluster_E[truthclustermaxlength] = {0};
  float m_truth_cluster_x[truthclustermaxlength] = {0};
  float m_truth_cluster_y[truthclustermaxlength] = {0};
  float m_truth_cluster_z[truthclustermaxlength] = {0};
  //looking at the code that generate the truth cluster, it seems that the time is always 0???
  float m_truth_cluster_t[truthclustermaxlength] = {0};
  int m_truth_cluster_detid[truthclustermaxlength] = {0};
  unsigned int m_truth_cluster_id[truthclustermaxlength] = {0};
  //the track id for the truth particle that generated the truth cluster
  int m_truth_cluster_trparticle_track_id[truthclustermaxlength] = {0};
  int m_nTruthClusters = 0;

  static const int truthtrackg4hitmaxlength = 1E7;
  float m_track_g4hit_x[truthtrackg4hitmaxlength] = {0};
  float m_track_g4hit_y[truthtrackg4hitmaxlength] = {0};
  float m_track_g4hit_z[truthtrackg4hitmaxlength] = {0};
  float m_track_g4hit_t[truthtrackg4hitmaxlength] = {0};
  float m_track_g4hit_E[truthtrackg4hitmaxlength] = {0};
  int m_track_g4hit_trparticle_track_id[truthtrackg4hitmaxlength] = {0};
  ULong64_t m_track_g4hit_id[truthtrackg4hitmaxlength] = {0};
  int m_nTrackG4Hits = 0;

  //sim stuff
  static const int ptruthmaxlength = 1E7;
  int m_particle_pid[ptruthmaxlength] = {0};
  float m_particle_energy[ptruthmaxlength] = {0};
  float m_particle_px[ptruthmaxlength] = {0};
  float m_particle_py[ptruthmaxlength] = {0};
  float m_particle_pz[ptruthmaxlength] = {0};
  float m_particle_vtx_x[ptruthmaxlength] = {0};
  float m_particle_vtx_y[ptruthmaxlength] = {0};
  float m_particle_vtx_z[ptruthmaxlength] = {0};
  int m_particle_track_id[ptruthmaxlength] = {0};
  int m_particle_primary_id[ptruthmaxlength] = {0};
  int m_particle_parent_id[ptruthmaxlength] = {0};
  int m_particle_is_pythia_primary[ptruthmaxlength] = {0};
  int m_particle_is_embedded[ptruthmaxlength] = {0};
  int m_particle_charge[ptruthmaxlength] = {0};
  int m_nParticles = 0;

   enum detid{
    mvtxId = 0,
    inttId = 1,
    tpcId = 2,
    tpotId = 3,
    cemcId = 4,
    ihcalId = 5,
    ohcalId = 6,
    epdId = 7,
    mbdId = 8
  };

  std::vector<std::pair<detid, TrkrDefs::TrkrId>> trkrlist = {
        {mvtxId, TrkrDefs::mvtxId},
        {inttId, TrkrDefs::inttId},
        {tpcId, TrkrDefs::tpcId},
        {tpotId, TrkrDefs::micromegasId}
    };

  SvtxEvalStack *m_svtxEvalStack = nullptr;
 
 
  TDatabasePDG *_pdg = nullptr;
  std::string Outfile;
  TFile *out;

};

#endif // calotrkana_H
