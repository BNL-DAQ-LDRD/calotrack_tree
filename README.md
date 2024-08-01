
# sPHENIX Simulation Result TTree Maker

This module creates a TTree containing the simulated detector response for all sPHENIX detectors (except for ZDC) and the primary truth particle information. 

## Detector Response Branches

The branches for the simulated detector responses combined for all detectors are as follows:

```cpp
T->Branch("nHits", &m_nHits, "nHits/I"); // Number of total detector hits
T->Branch("Hit_E", &m_Hit_E, "Hit_E[nHits]/F"); // Simulated energy
T->Branch("Hit_x", &m_Hit_x, "Hit_x[nHits]/F"); // Simulated spatial location [cm]
T->Branch("Hit_y", &m_Hit_y, "Hit_y[nHits]/F");
T->Branch("Hit_z", &m_Hit_z, "Hit_z[nHits]/F");
T->Branch("Hit_t", &m_Hit_t, "Hit_t[nHits]/F"); // Time of the hit [ns]
T->Branch("Hit_detid", &m_Hit_detid, "Hit_detid[nHits]/I"); // Detector ID
```

Hits for the tracking system (MVTX, INTT, TPC, TPOT) are represented by the `TrkrCluster` object, and hits for calorimeters (CEMC, Inner HCal, Outer HCal, sEPD) are represented by the `TowerInfo` object. Hits for the MBD are represented by the `MbdPmtHit` object.

Since there is no hit timing information for the tracking system, `Hit_t` is set to 0 for all tracking hits.

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
T->Branch("particle_energy", &m_particle_energy, "particle_energy[nParticles]/F"); // Particle energy[GeV]
T->Branch("particle_px", &m_particle_px, "particle_px[nParticles]/F"); // Particle momentum in x-direction [GeV/c]
T->Branch("particle_py", &m_particle_py, "particle_py[nParticles]/F"); // Particle momentum in y-direction
T->Branch("particle_pz", &m_particle_pz, "particle_pz[nParticles]/F"); // Particle momentum in z-direction
T->Branch("particle_vtx_x", &m_particle_vtx_x, "particle_vtx_x[nParticles]/F"); // Particle production vertex x-coordinate [cm]
T->Branch("particle_vtx_y", &m_particle_vtx_y, "particle_vtx_y[nParticles]/F"); // Particle production vertex y-coordinate
T->Branch("particle_vtx_z", &m_particle_vtx_z, "particle_vtx_z[nParticles]/F"); // Particle production vertex z-coordinate
```

