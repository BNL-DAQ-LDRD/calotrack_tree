

#include <G4_ActsGeom.C>


#include <calowaveformsim/CaloWaveformSim.h>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>

#include <fun4allraw/Fun4AllPrdfInputManager.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <globalvertex/GlobalVertexReco.h>

#include <phool/recoConsts.h>

#include <calotrkana/calotrkana.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libCaloWaveformSim.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libcalotrkana.so)


void Fun4All_run_dst(
    const int nEvents = 100,
    const string &inputFile0 = "g4hits.list",
    //const string &inputFile1 = "dst_calo_cluster.list",
    const string &inputFile1 = "dst_calo_waveform.list",
    const string &inputFile2 = "dst_trkr_cluster.list",
    const string &inputFile3 = "dst_mbd_epd.list",    
    const string &outputFile = "output_sim.root",

    const string &outputDSTFile = "DST_CALO_WAVEFORM_pp-0000000011-00000.root",
    const string &outDSTdir = ".",
    const string &cdbtag = "ProdA_2024")
{

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);
  recoConsts *rc = recoConsts::instance();

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", cdbtag);
  rc->set_uint64Flag("TIMESTAMP", 15);
  CDBInterface::instance()->Verbosity(1);

  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);

  Fun4AllInputManager *g4In = new Fun4AllDstInputManager("g4truth");
  g4In->AddListFile(inputFile0,1);
  se->registerInputManager(g4In);

  Fun4AllInputManager *CaloIn = new Fun4AllDstInputManager("calo");
  CaloIn->AddListFile(inputFile1,1);
  se->registerInputManager(CaloIn);

  Fun4AllInputManager *TrkrIn = new Fun4AllDstInputManager("trkr");
  TrkrIn->AddListFile(inputFile2,1);
  se->registerInputManager(TrkrIn);

  Fun4AllInputManager *MBDEPDIn = new Fun4AllDstInputManager("mbdepd");
  MBDEPDIn->AddListFile(inputFile3,1);
  se->registerInputManager(MBDEPDIn);

  Fun4AllInputManager *intrue2 = new Fun4AllRunNodeInputManager("DST_GEO");
  std::string geoLocation = CDBInterface::instance()->getUrl("calo_geo");
  intrue2->AddFile(geoLocation);
  se->registerInputManager(intrue2);

  ACTSGEOM::ActsGeomInit();

  calotrkana *caloana24 = new calotrkana("calotrkana","testout.root");
  se->registerSubsystem(caloana24);

  se->run(nEvents);
  CDBInterface::instance()->Print();  // print used DB files
  se->End();
  cout << "JOB COMPLETE." <<endl;
  se->PrintTimer();
  gSystem->Exit(0);
}
