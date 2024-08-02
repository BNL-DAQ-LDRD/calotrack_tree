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
  static const int recomaxlength = 1E6;
  
  float m_Hit_E[recomaxlength] = {0};
  float m_Hit_x[recomaxlength] = {0};
  float m_Hit_y[recomaxlength] = {0};
  float m_Hit_z[recomaxlength] = {0};
  float m_Hit_t[recomaxlength] = {0};
  int m_Hit_detid[recomaxlength] = {0};
  int m_nHits = 0;
  //sim stuff
  static const int ptruthmaxlength = 1E5;
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
 
 

  std::string Outfile;
  TFile *out;

};

#endif // calotrkana_H
