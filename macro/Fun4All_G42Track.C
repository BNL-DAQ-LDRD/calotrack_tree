#ifndef MACRO_FUN4ALLG4SPHENIX_C
#define MACRO_FUN4ALLG4SPHENIX_C

#include <GlobalVariables.C>

#include <DisplayOn.C>
#include <G4Setup_sPHENIX.C>
#include <G4_Mbd.C>
#include <G4_CaloTrigger.C>
#include <G4_Centrality.C>
#include <G4_DSTReader.C>
#include <G4_Global.C>
#include <G4_HIJetReco.C>
#include <G4_Input.C>
#include <G4_Jets.C>
#include <G4_KFParticle.C>
#include <G4_ParticleFlow.C>
#include <G4_Production.C>
#include <G4_TopoClusterReco.C>

#include <Trkr_RecoInit.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>
#include <Trkr_Reco.C>
#include <Trkr_Eval.C>
#include <Trkr_QA.C>

#include <Trkr_Diagnostics.C>
#include <G4_User.C>
#include <QA.C>

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <calotrkana/calotrkana.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)

// For HepMC Hijing
// try inputFile = /sphenix/sim/sim01/sphnxpro/sHijing_HepMC/sHijing_0-12fm.dat

int Fun4All_G42Track(
    const int nEvents = 1,
    const string &inputFile = "/sphenix/user/jdosbo/fm4npp/macros/detectors/sPHENIX/g4hits.list",
    const string &outputFile = "G4sPHENIX.root",
    const string &embed_input_file = "https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/sPHENIX_G4Hits_sHijing_9-11fm_00000_00010.root",
    const int skip = 0,
    const string &outdir = ".")
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);

  //Opt to print all random seed used for debugging reproducibility. Comment out to reduce stdout prints.
  PHRandomSeed::Verbosity(1);
  CDBInterface::instance()->Verbosity(1);
  CDB::global_tag="MDC2_ana.509";
  
  // just if we set some flags somewhere in this macro
  recoConsts *rc = recoConsts::instance();
  // By default every random number generator uses
  // PHRandomSeed() which reads /dev/urandom to get its seed
  // if the RANDOMSEED flag is set its value is taken as seed
  // You can either set this to a random value using PHRandomSeed()
  // which will make all seeds identical (not sure what the point of
  // this would be:
  //  rc->set_IntFlag("RANDOMSEED",PHRandomSeed());
  // or set it to a fixed value so you can debug your code
  //  rc->set_IntFlag("RANDOMSEED", 12345);

  ifstream infile(inputFile);
  string line;
  getline(infile, line);
  rc->set_StringFlag("Sim_File_Name", line);

  Input::READHITS = true;
  INPUTREADHITS::listfile[0] = inputFile;

  //-----------------
  // Initialize the selected Input/Event generation
  //-----------------
  // This creates the input generator(s)
  InputInit();

  // register all input generators with Fun4All
  InputRegister();

  if (! Input::READHITS)
  {
    rc->set_IntFlag("RUNNUMBER",1);

    SyncReco *sync = new SyncReco();
    se->registerSubsystem(sync);

    HeadReco *head = new HeadReco();
    se->registerSubsystem(head);
  }
// Flag Handler is always needed to read flags from input (if used)
// and update our rc flags with them. At the end it saves all flags
// again on the DST in the Flags node under the RUN node
  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);
  
  Enable::TRACKING_VERBOSITY = 0;
  
  Enable::PIPE = true;
  Enable::PIPE_ABSORBER = true;

  // central tracking
  Enable::MVTX = true;
  Enable::MVTX_CELL = Enable::MVTX && true;
  Enable::MVTX_CLUSTER = Enable::MVTX_CELL && true;

  Enable::INTT = true;
//  Enable::INTT_ABSORBER = true; // enables layerwise support structure readout
//  Enable::INTT_SUPPORT = true; // enable global support structure readout
  Enable::INTT_CELL = Enable::INTT && true;
  Enable::INTT_CLUSTER = Enable::INTT_CELL && true;

  Enable::TPC = true;
  Enable::TPC_ABSORBER = true;
  Enable::TPC_CELL = Enable::TPC && true;
  Enable::TPC_CLUSTER = Enable::TPC_CELL && true;

  Enable::MICROMEGAS = true;
  Enable::MICROMEGAS_CELL = Enable::MICROMEGAS && true;
  Enable::MICROMEGAS_CLUSTER = Enable::MICROMEGAS_CELL && true;

  Enable::TRACKING_TRACK = (Enable::MICROMEGAS_CLUSTER && Enable::TPC_CLUSTER && Enable::INTT_CLUSTER && Enable::MVTX_CLUSTER) && true;
  Enable::GLOBAL_RECO = (Enable::MBDFAKE || Enable::MBDRECO || Enable::TRACKING_TRACK) && true;
  Enable::TRACKING_EVAL = Enable::TRACKING_TRACK && Enable::GLOBAL_RECO && true;

  Enable::MAGNET = true;
  Enable::MAGNET_ABSORBER = true;


  Enable::CDB = true;
  // global tag
  rc->set_StringFlag("CDB_GLOBALTAG",CDB::global_tag);
  // 64 bit timestamp
  rc->set_uint64Flag("TIMESTAMP",CDB::timestamp);

  G4Init();

  //---------------------
  // GEANT4 Detector description
  //---------------------
  if (!Input::READHITS)
  {
    G4Setup();
  }

  //------------------
  // Detector Division
  //------------------

  if (Enable::MVTX_CELL) Mvtx_Cells();
  if (Enable::INTT_CELL) Intt_Cells();
  if (Enable::TPC_CELL) TPC_Cells();
  if (Enable::MICROMEGAS_CELL) Micromegas_Cells();

  //--------------
  // SVTX tracking
  //--------------
  if(Enable::TRACKING_TRACK)
    {
      TrackingInit();
    }
  if (Enable::MVTX_CLUSTER) Mvtx_Clustering();
  if (Enable::INTT_CLUSTER) Intt_Clustering();
  if (Enable::TPC_CLUSTER)
    {
      if(G4TPC::ENABLE_DIRECT_LASER_HITS || G4TPC::ENABLE_CENTRAL_MEMBRANE_HITS)
	{
	  TPC_LaserClustering();
	}
      else
	{
	  TPC_Clustering();
	}
    }
  if (Enable::MICROMEGAS_CLUSTER) Micromegas_Clustering();

  if (Enable::TRACKING_TRACK)
  {
    Tracking_Reco_TrackSeed();
    
  auto converter = new TrackSeedTrackMapConverter;
  // Default set to full SvtxTrackSeeds. Can be set to
  // SiliconTrackSeedContainer or TpcTrackSeedContainer
  converter->setTrackSeedName("TpcTrackSeedContainer");
  converter->setFieldMap(G4MAGNET::magfield_tracking);
  converter->Verbosity(0);
  se->registerSubsystem(converter);

  se->registerSubsystem(new PHSimpleVertexFinder);
  }


  //-----------------
  // Global Vertexing
  //-----------------

  if (Enable::GLOBAL_RECO && Enable::GLOBAL_FASTSIM)
  {
    cout << "You can only enable Enable::GLOBAL_RECO or Enable::GLOBAL_FASTSIM, not both" << endl;
    gSystem->Exit(1);
  }
  if (Enable::GLOBAL_RECO)
  {
    Global_Reco();
  }
  else if (Enable::GLOBAL_FASTSIM)
  {
    Global_FastSim();
  }

  //----------------------
  // Simulation evaluation
  //----------------------
  string outputroot = outputFile;
  string remove_this = ".root";
  size_t pos = outputroot.find(remove_this);
  if (pos != string::npos)
  {
    outputroot.erase(pos, remove_this.length());
  }

  if (Enable::TRACKING_EVAL) Tracking_Eval(outputroot + "_g4svtx_eval.root");

  //add calotrkana code here

  

  //--------------
  // Set up Input Managers
  //--------------

  InputManagers();


  // if we use a negative number of events we go back to the command line here
  if (nEvents < 0)
  {
    return 0;
  }
  // if we run the particle generator and use 0 it'll run forever
  // for embedding it runs forever if the repeat flag is set
  if (nEvents == 0 && !Input::HEPMC && !Input::READHITS && INPUTEMBED::REPEAT)
  {
    cout << "using 0 for number of events is a bad idea when using particle generators" << endl;
    cout << "it will run forever, so I just return without running anything" << endl;
    return 0;
  }

  calotrkana *caloana24 = new calotrkana("calotrkana", "calotrkana.root");
  caloana24->set_tracking_only(true);
  // caloana24->set_npart_range(0, 32);
  se->registerSubsystem(caloana24);

  se->skip(skip);
  se->run(nEvents);
  //  se->PrintTimer();

  //-----
  // Exit
  //-----

  CDBInterface::instance()->Print(); // print used DB files
  se->End();
  std::cout << "All done" << std::endl;
  delete se;


  gSystem->Exit(0);
  return 0;
}
#endif
