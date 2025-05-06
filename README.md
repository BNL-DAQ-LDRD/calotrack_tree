# sPHENIX Simulation Result TTree Maker

This module creates a TTree containing the simulated detector response for all sPHENIX detectors (except for ZDC) and the primary truth particle information.

## cheat sheet

Install:
ref: https://github.com/sPHENIX-Collaboration/tutorials/blob/master/CreateSubsysRecoModule/source/create_me.sh

Run:
```bash
root -l -q Fun4All_run_dst.C >& log
root -l -q Fun4All_G4_sPHENIX.C\(1\) >& log
```


## Detector Response Branches

The branches for the simulated detector responses for calorimeters and MBD detectors are as follows:

```cpp
T->Branch("nHits", &m_nHits, "nHits/I"); // Number of total detector hits
T->Branch("Hit_E", &m_Hit_E, "Hit_E[nHits]/F"); // Simulated energy
T->Branch("Hit_x", &m_Hit_x, "Hit_x[nHits]/F"); // Simulated spatial location [cm]
T->Branch("Hit_y", &m_Hit_y, "Hit_y[nHits]/F");
T->Branch("Hit_z", &m_Hit_z, "Hit_z[nHits]/F");
T->Branch("Hit_t", &m_Hit_t, "Hit_t[nHits]/F"); // Time of the hit [ns]
T->Branch("Hit_detid", &m_Hit_detid, "Hit_detid[nHits]/I"); // Detector ID
```

Hits for the calorimeters (CEMC, Inner HCal, Outer HCal, sEPD) are represented by the `TowerInfo` object. Hits for the MBD are represented by the `MbdPmtHit` object.

### Detector ID Definition

The detector ID is defined as follows:

```cpp
enum detid {
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
```

## Primary Truth Particle Information Branches

The branches for the primary truth particle information are as follows:

```cpp
T->Branch("nParticles", &m_nParticles, "nParticles/I"); // Number of primary particles
T->Branch("particle_pid", &m_particle_pid, "particle_pid[nParticles]/I"); // Particle PID
T->Branch("particle_energy", &m_particle_energy, "particle_energy[nParticles]/F"); // Particle energy [GeV]
T->Branch("particle_px", &m_particle_px, "particle_px[nParticles]/F"); // Particle momentum in x-direction [GeV/c]
T->Branch("particle_py", &m_particle_py, "particle_py[nParticles]/F"); // Particle momentum in y-direction
T->Branch("particle_pz", &m_particle_pz, "particle_pz[nParticles]/F"); // Particle momentum in z-direction
T->Branch("particle_vtx_x", &m_particle_vtx_x, "particle_vtx_x[nParticles]/F"); // Particle production vertex x-coordinate [cm]
T->Branch("particle_vtx_y", &m_particle_vtx_y, "particle_vtx_y[nParticles]/F"); // Particle production vertex y-coordinate
T->Branch("particle_vtx_z", &m_particle_vtx_z, "particle_vtx_z[nParticles]/F"); // Particle production vertex z-coordinate
T->Branch("particle_track_id", &m_particle_track_id, "particle_track_id[nParticles]/I"); // the unique track id for the particle
T->Branch("particle_primary_id", &m_particle_primary_id, "particle_primary_id[nParticles]/I");// the track id for the primary particle that result in current particle
T->Branch("particle_parent_id", &m_particle_parent_id, "particle_parent_id[nParticles]/I"); // the parent particle track id
T->Branch("particle_is_pythia_primary", &m_particle_is_pythia_primary, "particle_is_pythia_primary[nParticles]/I"); // if the particle is primary particle
T->Branch("particle_charge", &m_particle_charge, "particle_charge[nParticles]/I"); //the particle electric charge
```

## Tracking System Information Branches

The branches for the tracking system information are as follows:

```cpp
T->Branch("nRecoClusters", &m_nRecoClusters, "nRecoClusters/I"); // Number of reconstructed clusters
T->Branch("reco_cluster_E", &m_reco_cluster_E, "reco_cluster_E[nRecoClusters]/F"); // Reconstructed cluster energy
T->Branch("reco_cluster_x", &m_reco_cluster_x, "reco_cluster_x[nRecoClusters]/F"); // Reconstructed cluster x-coordinate
T->Branch("reco_cluster_y", &m_reco_cluster_y, "reco_cluster_y[nRecoClusters]/F"); // Reconstructed cluster y-coordinate
T->Branch("reco_cluster_z", &m_reco_cluster_z, "reco_cluster_z[nRecoClusters]/F"); // Reconstructed cluster z-coordinate
T->Branch("reco_cluster_t", &m_reco_cluster_t, "reco_cluster_t[nRecoClusters]/F"); // Reconstructed cluster time
T->Branch("reco_cluster_detid", &m_reco_cluster_detid, "reco_cluster_detid[nRecoClusters]/I"); // Reconstructed cluster detector ID
T->Branch("reco_cluster_trcluster_id", &m_reco_cluster_trcluster_id, "reco_cluster_trcluster_id[nRecoClusters]/i"); // The best matching truth cluster id, 0 if not matched
T->Branch("reco_cluster_g4hit_id", &m_reco_cluster_g4hit_id, "reco_cluster_g4hit_id[nRecoClusters]/l"); // best matching G4Hit id, 0 if not matched (this id is a type unsigned long long)

T->Branch("nTruthClusters", &m_nTruthClusters, "nTruthClusters/I"); // Number of truth clusters
T->Branch("truth_cluster_E", &m_truth_cluster_E, "truth_cluster_E[nTruthClusters]/F"); // Truth cluster energy
T->Branch("truth_cluster_x", &m_truth_cluster_x, "truth_cluster_x[nTruthClusters]/F"); // Truth cluster x-coordinate
T->Branch("truth_cluster_y", &m_truth_cluster_y, "truth_cluster_y[nTruthClusters]/F"); // Truth cluster y-coordinate
T->Branch("truth_cluster_z", &m_truth_cluster_z, "truth_cluster_z[nTruthClusters]/F"); // Truth cluster z-coordinate
T->Branch("truth_cluster_t", &m_truth_cluster_t, "truth_cluster_t[nTruthClusters]/F"); // Truth cluster time
T->Branch("truth_cluster_detid", &m_truth_cluster_detid, "truth_cluster_detid[nTruthClusters]/I"); // Truth cluster detector ID
T->Branch("truth_cluster_id", &m_truth_cluster_id, "truth_cluster_id[nTruthClusters]/i"); //Unique truth cluster ID
T->Branch("truth_cluster_trparticle_track_id", &m_truth_cluster_trparticle_track_id, "truth_cluster_trparticle_track_id[nTruthClusters]/I"); // Truth particle track ID associated with the truth cluster
```
Detector id for traking is defined [here](https://github.com/sPHENIX-Collaboration/coresoftware/blob/master/offline/packages/trackbase/TrkrDefs.h#L51-L57):
```
enum TrkrId
  {
    mvtxId = 0,
    inttId = 1,
    tpcId = 2,
    micromegasId = 3
  };
```
8/24 update
Now we also include the matching G4Hit information
```cpp
T->Branch("nTrackG4Hits", &m_nTrackG4Hits, "nTrackG4Hits/I");
T->Branch("track_g4hit_x", &m_track_g4hit_x, "track_g4hit_x[nTrackG4Hits]/F");
T->Branch("track_g4hit_y", &m_track_g4hit_y, "track_g4hit_y[nTrackG4Hits]/F");
T->Branch("track_g4hit_z", &m_track_g4hit_z, "track_g4hit_z[nTrackG4Hits]/F");
T->Branch("track_g4hit_t", &m_track_g4hit_t, "track_g4hit_t[nTrackG4Hits]/F");
T->Branch("track_g4hit_E", &m_track_g4hit_E, "track_g4hit_E[nTrackG4Hits]/F");// G4Hit total energy deposition 
T->Branch("track_g4hit_trparticle_track_id", &m_track_g4hit_trparticle_track_id, "track_g4hit_trparticle_track_id[nTrackG4Hits]/I");// the track id of the G4particle that creates this hit
T->Branch("track_g4hit_id", &m_track_g4hit_id, "track_g4hit_id[nTrackG4Hits]/l"); //unique G4Hit id
```

Some truth matching is done for the tracking information thanks to the [g4eval lib](https://github.com/sPHENIX-Collaboration/coresoftware/tree/master/simulation/g4simulation/g4eval), each reco cluster(object equivalent to what we get in real data) is matched to a truth cluster(a set of G4Hits) by the unique truth cluster id(0 if no matching is found). Each truth cluster is linked to a truth particle(int the `particle_` branches) by the particle track id.


## Some Notes:
Currently the truth matching is not working correctly with the existing DST files due to the fact that some nodes are produced under different ana build, and I suspect the different in traking detector geom is messing up the truth cluster matching. For example, while using our run 15 pythia8 pp production, none of the reco cluster in the MVTX and INTT can find any matched truth cluster. This problem is gone after me running the simulation with a consistant build from the GEANT4 step.(something we need to keep in mind while using the produced files). I also noticed that the truth association is turned off for our current production, maybe we can request for the new sPHENIX MB simulation production(with the new tune) to have different passes produced under the same ana build and (maybe) have the TPC truth association on during the production. 
