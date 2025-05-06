import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from util import *
from sklearn.metrics import adjusted_rand_score
from util import *

aris = []

bin_edges = np.linspace(0, 5, 26)  # Define bin edges
total_matched_hist = np.zeros(len(bin_edges) - 1, dtype=int)
total_all_hist = np.zeros(len(bin_edges) - 1, dtype=int)

# Loop to read and process files from event_0 to event_10
for ievent in range(2,3):
    print (f"Processing event {ievent}...")
    with pd.HDFStore(f'2025-05-05-pp-full-seeding-100event/data_event_{ievent}.h5', mode='r') as store:
        clusters = store['clusters']      # Saved as fixed format
        seeds = store['seeds']            # Saved as fixed format
        particles = store['particles']    # Saved as fixed format

        print(f"Number of particles: {len(particles)}")
        particles = particles[particles['eta'].apply(abs) < 1.1]
        print(f"|eta| < 1.1: {len(particles)}")
        particles = particles[particles['vz'].apply(abs) < 10]
        print(f"|vz| < 10: {len(particles)}")
        particles = particles[particles['ptid']>0]
        print(f"ptid > 0: {len(particles)}")

        ncommon = 30  # change to your desired threshold
        matched_pt = match_particles_to_seeds_optimized(particles, seeds, ncommon)
        all_pt = particles['pt'].tolist()

        matched_hist, _ = np.histogram(matched_pt, bins=bin_edges)
        all_hist, _ = np.histogram(all_pt, bins=bin_edges)
        total_matched_hist += matched_hist
        total_all_hist += all_hist
        print(f"Event {ievent}: All {len(all_pt)} matched: {len(matched_pt)}\n")

        # Calculate ARI
        cid_to_index = {cid: index for index, cid in enumerate(clusters['cid'])}
        groupids_particle = get_group_ids(particles, clusters.shape[0], cid_to_index)
        print(f"groupids_particle {groupids_particle}")
        groupids_seed = get_group_ids(seeds, clusters.shape[0], cid_to_index)
        print(f"groupids_seed {groupids_seed}")
        ari = adjusted_rand_score(groupids_seed, groupids_particle)
        aris.append(ari)
        print(f"Event {ievent} ARI: {ari}")

print("Total Matched histogram:", total_matched_hist)
print("Total All histogram:", total_all_hist)
print("average ARI:", np.mean(aris))

# Plotting the histograms
plot_eff((total_all_hist, bin_edges), (total_matched_hist, bin_edges))