

#include <G4_ActsGeom.C>

#include <GlobalVariables.C>

#include <G4_Production.C>
#include <G4_Centrality.C>
#include <Trkr_RecoInit.C>
#include <Trkr_Clustering.C>

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

#include <centrality/CentralityReco.h>

//tracking
#include <trackreco/PHCASeeding.h>

#include <trackingdiagnostics/TrackSeedTrackMapConverter.h>

#include <calotrkana/calotrkana.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libCaloWaveformSim.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libcalotrkana.so)
R__LOAD_LIBRARY(libcentrality_io.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)
R__LOAD_LIBRARY(libtrack_reco.so)

void Fun4All_run_dst(
    const int nEvents = 10,
    const string &inputFile0 = "g4hits.list",
    // const string &inputFile1 = "dst_global.list",
    const string &inputFile1 = "dst_calo_waveform.list",
    const string &inputFile2 = "dst_global.list",
    const string &inputFile3 = "dst_mbd_epd.list",
    const string &inputFile4 = "dst_tracks.list",
    const string &inputFile5 = "dst_trkr_hit.list",
    const string &inputFile6 = "dst_truth.list",
    const string &outputFile = "output_sim.root",

    const string &outputDSTFile = "DST_CALO_WAVEFORM_pp-0000000011-00000.root",
    const string &outDSTdir = ".",
    const string &cdbtag = "MDC2_ana.435")
{

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(2);
  recoConsts *rc = recoConsts::instance();

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", cdbtag);
  rc->set_uint64Flag("TIMESTAMP", 19);
  //get the first line of inputFile0
  ifstream infile(inputFile0);
  string line;
  getline(infile, line);
  rc->set_StringFlag("Sim_File_Name", line);
  CDBInterface::instance()->Verbosity(1);

  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);

  Fun4AllInputManager *CaloIn = new Fun4AllDstInputManager("calo");
  CaloIn->AddListFile(inputFile1, 1);
  se->registerInputManager(CaloIn);

  Fun4AllInputManager *GlobIn = new Fun4AllDstInputManager("global");
  GlobIn->AddListFile(inputFile2, 1);
  se->registerInputManager(GlobIn);

  Fun4AllInputManager *MBDEPDIn = new Fun4AllDstInputManager("mbdepd");
  MBDEPDIn->AddListFile(inputFile3, 1);
  se->registerInputManager(MBDEPDIn);
  /*
  //for the eval code
  Fun4AllInputManager *TrackIn = new Fun4AllDstInputManager("track");
  TrackIn->AddListFile(inputFile4,1);
  se->registerInputManager(TrackIn);
*/
  Fun4AllInputManager *TrkrHitIn = new Fun4AllDstInputManager("trkrhit");
  TrkrHitIn->AddListFile(inputFile5, 1);
  se->registerInputManager(TrkrHitIn);

  Fun4AllInputManager *g4In = new Fun4AllDstInputManager("g4truth");
  g4In->AddListFile(inputFile0, 1);
  se->registerInputManager(g4In);

  Fun4AllInputManager *TruthIn = new Fun4AllDstInputManager("truth");
  TruthIn->AddListFile(inputFile6, 1);
  se->registerInputManager(TruthIn);

  // need to redo TPC digitization for the truth association

  Fun4AllInputManager *intrue2 = new Fun4AllRunNodeInputManager("DST_GEO");
  std::string geoLocation = CDBInterface::instance()->getUrl("calo_geo");
  intrue2->AddFile(geoLocation);
  se->registerInputManager(intrue2);

  // central tracking
  Enable::MVTX = true;
  Enable::INTT = true;
  Enable::TPC = true;
  Enable::MICROMEGAS = true;

  // TPC
  G4TPC::ENABLE_STATIC_DISTORTIONS = false;
  G4TPC::ENABLE_TIME_ORDERED_DISTORTIONS = false;
  // G4TPC::ENABLE_CORRECTIONS = false;
  G4TPC::DO_HIT_ASSOCIATION = true;

  // tracking configuration
  G4TRACKING::use_full_truth_track_seeding = false;

  // do not initialize magnetic field in ACTS
  G4TRACKING::init_acts_magfield = false;

  MagnetFieldInit();
  TrackingInit();

  // clustering
  Mvtx_Clustering();
  Intt_Clustering();
  TPC_Clustering();
  Micromegas_Clustering();

  // TPC seeder
  {
    int verbosity = std::max(Enable::VERBOSITY, Enable::TRACKING_VERBOSITY);
    auto seeder = new PHCASeeding("PHCASeeding");
    double fieldstrength = std::numeric_limits<double>::quiet_NaN(); // set by isConstantField if constant
    bool ConstField = isConstantField(G4MAGNET::magfield_tracking, fieldstrength);
    if (ConstField)
    {
      seeder->useConstBField(true);
      seeder->constBField(fieldstrength);
    }
    else
    {
      seeder->set_field_dir(-1 * G4MAGNET::magfield_rescale);
      seeder->useConstBField(false);
      seeder->magFieldFile(G4MAGNET::magfield_tracking); // to get charge sign right
    }
    seeder->Verbosity(verbosity);
    seeder->SetLayerRange(7, 55);
    seeder->SetSearchWindow(1.5, 0.05); // (z width, phi width)
    seeder->SetMinHitsPerCluster(0);
    seeder->SetMinClustersPerTrack(3);
    seeder->useFixedClusterError(true);

    seeder->set_pp_mode(true);
    se->registerSubsystem(seeder);
  }

  {
    int verbosity = std::max(Enable::VERBOSITY, Enable::TRACKING_VERBOSITY);

    auto converter = new TrackSeedTrackMapConverter;
    // Default set to full SvtxTrackSeeds. Can be set to
    // SiliconTrackSeedContainer or TpcTrackSeedContainer
    converter->setTrackSeedName("TpcTrackSeedContainer");
    converter->setFieldMap(G4MAGNET::magfield_tracking);
    converter->Verbosity(verbosity);
    se->registerSubsystem(converter);
  }
  //this is in TrackingInit();
  //ACTSGEOM::ActsGeomInit();

  // Centrality();

  calotrkana *caloana24 = new calotrkana("calotrkana", "testout.root");
  caloana24->set_npart_range(0, 32);
  se->registerSubsystem(caloana24);

  se->run(nEvents);
  CDBInterface::instance()->Print(); // print used DB files
  se->End();
  cout << "JOB COMPLETE." << endl;
  se->PrintTimer();
  gSystem->Exit(0);
}
