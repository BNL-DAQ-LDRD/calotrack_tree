import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from util import *
from sklearn.metrics import adjusted_rand_score
from util import *


bin_edges = np.linspace(0, 5, 26)  # Define bin edges
total_matched_hist = np.zeros(len(bin_edges) - 1, dtype=int)
total_all_hist = np.zeros(len(bin_edges) - 1, dtype=int)
binned_aris = [ [] for _ in range(len(bin_edges) - 1)]

# Loop to read and process files from event_0 to event_10
for ievent in range(0, 100):
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
        cid_to_pt = cid_to_pt_mapping(particles)
        groupids_particle = get_group_ids(particles, clusters.shape[0], cid_to_index, cid_to_pt, bin_edges)
        print(f"groupids_particle {groupids_particle}")
        groupids_seed = get_group_ids(seeds, clusters.shape[0], cid_to_index, cid_to_pt, bin_edges)
        print(f"groupids_seed {groupids_seed}")
        for ibin in range(len(bin_edges) - 1):
            ari = adjusted_rand_score(groupids_particle[ibin], groupids_seed[ibin])
            binned_aris[ibin].append(ari)
            print(f"ievent {ievent}, bin {ibin}: ARI = {ari:.4f}")

print("Total Matched histogram:", total_matched_hist)
print("Total All histogram:", total_all_hist)
# Calculate average ARI for each bin
avg_aris = [np.mean(bin_aris) if bin_aris else 0 for bin_aris in binned_aris]
std_aris = [np.std(bin_aris) if len(bin_aris) > 1 else 0 for bin_aris in binned_aris]

# Get bin centers for plotting
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# Plot average ARI by pT bin
plt.figure(figsize=(10, 6))
plt.errorbar(bin_centers, avg_aris, yerr=std_aris, marker='o', linestyle='-', capsize=3)
plt.xlabel('pT (GeV/c)')
plt.ylabel('Average Adjusted Rand Index')
plt.title('Clustering Quality (ARI) vs. Particle pT')
plt.grid(True, alpha=0.3)
plt.savefig('ari_vs_pt.png')
plt.show()

# Plotting the histograms
plot_eff((total_all_hist, bin_edges), (total_matched_hist, bin_edges))