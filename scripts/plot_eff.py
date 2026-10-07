import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from util import *
from sklearn.metrics import adjusted_rand_score
from util import *
from tqdm import tqdm
from matplotlib.colors import LogNorm

bin_edges = np.linspace(0, 5, 26)  # Define bin edges
total_matched_hist = np.zeros(len(bin_edges) - 1, dtype=int)
total_all_hist = np.zeros(len(bin_edges) - 1, dtype=int)
total_matched_seed_hist = np.zeros(len(bin_edges) -1, dtype=int)
total_all_seed_hist = np.zeros(len(bin_edges) -1, dtype=int)
all_truth_pt = []
all_reco_pt = []

all_reco_ntpc = []

binned_aris = [ [] for _ in range(len(bin_edges) - 1)]

integrated_pt_bins = np.array([0.5, 1.0, 1000.0])  # Example bin edges for pT
integrated_pt_matched_hist = np.zeros(len(integrated_pt_bins) - 1, dtype=int)
integrated_pt_all_hist = np.zeros(len(integrated_pt_bins) - 1, dtype=int)

dummy_pt_bin_edges = np.linspace(0, 1000, 2)  # a dummy bin to include all particles
unbinned_aris = [ [] for _ in range(len(dummy_pt_bin_edges) - 1)]
# Initialize tqdm for the event processing loop
event_range = range(0, 9000)
progress_bar = tqdm(event_range, desc="Processing events", unit="event")
# Loop to read and process files from event_0 to event_10
for ievent in progress_bar:
    # print (f"Processing event {ievent}...")
    # fname = f'2025-05-13-pp-1k-cut/data_event_{ievent}.h5'
    # fname = f'2025-05-20-pp-1k-nocut/data_event_{ievent}.h5'
    fname = f'2025-05-20-pp-10k-mincut/data_event_{ievent}.h5'
    with pd.HDFStore(fname, mode='r') as store:
        clusters = store['clusters']      # Saved as fixed format
        seeds = store['seeds']            # Saved as fixed format
        particles = store['particles']    # Saved as fixed format

        # print(f"Number of particles: {len(particles)}")
        particles = particles[particles['cids'].apply(len) > 5]
        # print(f"|nclus| > 20: {len(particles)}")
        particles = particles[particles['eta'].apply(abs) < 1.1]
        # print(f"|eta| < 1.1: {len(particles)}")
        # particles = particles[particles['vz'].apply(abs) < 10]
        # print(f"|vz| < 10: {len(particles)}")
        # particles = particles[particles['ptid']>0]
        # print(f"ptid > 0: {len(particles)}")

        ncommon = 24  # change to your desired threshold
        matched_pt = match_particles_to_seeds_optimized(particles, seeds, ncommon)
       
        all_pt = particles['pt'].tolist()
        
        all_reco_ntpc.extend((seeds['cids'].apply(len)).tolist())
        
        matched_hist, _ = np.histogram(matched_pt, bins=bin_edges)
        all_hist, _ = np.histogram(all_pt, bins=bin_edges)
        total_matched_hist += matched_hist
        total_all_hist += all_hist
        # print(f"Event {ievent}: All {len(all_pt)} matched: {len(matched_pt)}\n")

        matched_hist, _ = np.histogram(matched_pt, bins=integrated_pt_bins)
        all_hist, _ = np.histogram(all_pt, bins=integrated_pt_bins)
        integrated_pt_matched_hist += matched_hist
        integrated_pt_all_hist += all_hist
        # print(f"Event {ievent}: All {len(integrated_pt_all_hist)} matched: {len(integrated_pt_matched_hist)}\n")

        # Calculate ARI
        cid_to_index = {cid: index for index, cid in enumerate(clusters['cid'])}
        cid_to_pt = cid_to_pt_mapping(particles)
        groupids_particle = get_group_ids(particles, clusters.shape[0], cid_to_index, cid_to_pt, bin_edges)
        # print(f"groupids_particle {groupids_particle}")
        groupids_seed = get_group_ids(seeds, clusters.shape[0], cid_to_index, cid_to_pt, bin_edges)
        # print(f"groupids_seed {groupids_seed}")
        for ibin in range(len(bin_edges) - 1):
            ari = adjusted_rand_score(groupids_particle[ibin], groupids_seed[ibin])
            binned_aris[ibin].append(ari)
            # print(f"ievent {ievent}, bin {ibin}: ARI = {ari:.4f}")
        
        groupids_particle = get_group_ids(particles, clusters.shape[0], cid_to_index, cid_to_pt, dummy_pt_bin_edges)
        groupids_seed = get_group_ids(seeds, clusters.shape[0], cid_to_index, cid_to_pt, dummy_pt_bin_edges)
        # print(f"groupids_particle {groupids_particle}")
        # print(f"groupids_seed {groupids_seed}")
        for ibin in range(len(dummy_pt_bin_edges) - 1):
            ari = adjusted_rand_score(groupids_particle[ibin], groupids_seed[ibin])
            unbinned_aris[ibin].append(ari)
            # print(f"ievent {ievent}, bin {ibin}: unbinned ARI = {ari:.4f}")

        
        # circle fit seeds to get purity
        if len(seeds) == 0:
            continue
        matched_seed_ids = match_seeds_to_particles_optimized(seeds, particles, ncommon)
        
        circle_fits = seeds.apply(fit_circle_row, axis=1)
        circle_fits['pt'] = 0.3 * 1.4 * (circle_fits['radius'] / 100)
        seeds['center_x'] = circle_fits['center_x']
        seeds['center_y'] = circle_fits['center_y']
        seeds['radius'] = circle_fits['radius']
        seeds['pt'] = circle_fits['pt']
        
        matched_seeds = seeds[seeds['matched'] == True]
        all_truth_pt.extend(matched_seeds['truthpt'].tolist())
        all_reco_pt.extend(matched_seeds['pt'].tolist())
        
        
        matched_pt_seeds = seeds[seeds['matched'] == True]['pt'].tolist()
        matched_hist_seed, _ = np.histogram(matched_pt_seeds, bins=bin_edges)
        total_matched_seed_hist += matched_hist_seed
        
        all_seed_pt = circle_fits['pt'].tolist()
        all_seed_hist, _ = np.histogram(all_seed_pt, bins=bin_edges)
        total_all_seed_hist += all_seed_hist


      
print("Total Matched histogram:", total_matched_hist)
print("Total All histogram:", total_all_hist)

print("Integrated pT Matched histogram:", integrated_pt_matched_hist)
print("Integrated pT All histogram:", integrated_pt_all_hist)

# Calculate average ARI for each bin
avg_aris = [np.mean(bin_aris) if bin_aris else 0 for bin_aris in binned_aris]
std_aris = [np.std(bin_aris) if len(bin_aris) > 1 else 0 for bin_aris in binned_aris]
print("avg_aris:", avg_aris)

# Calculate average and std for unbinned ARI
unbinned_avg_aris = [np.mean(bin_aris) if bin_aris else 0 for bin_aris in unbinned_aris]
unbinned_std_aris = [np.std(bin_aris) if len(bin_aris) > 1 else 0 for bin_aris in unbinned_aris]
print("unbinned_avg_aris", unbinned_avg_aris)


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

plt.figure(figsize=(10,6))
countsntpc, bin_edgesntpc, _ = plt.hist(all_reco_ntpc, bins=50, range=(0,50))
plt.xlabel('nTPC')
plt.show()
np.savez('ntpc_caseed_data.npz', counts=countsntpc, bin_edges=bin_edgesntpc)

bins = [np.linspace(0, 3, 51), np.linspace(0, 3, 51)]  # 50 bins from 0-3 GeV

plt.figure(figsize=(10, 8))
plt.hist2d(all_truth_pt, all_reco_pt, bins=bins, cmap='viridis', cmin=1, norm=LogNorm(vmin=0.1))
plt.xlabel('Truth p$_{T}$ [GeV/c]')
plt.ylabel('Reconstructed p$_{T}$ [GeV/c]')
plt.colorbar(label='')
plt.plot([0, max(all_truth_pt)], [0, max(all_truth_pt)], 'r--', lw=2)
plt.xlim(0, 3)
plt.ylim(0, 3)
plt.show()

# Plotting the histograms
plot_eff((total_all_hist, bin_edges), (total_matched_hist, bin_edges), "caseed_eff.npz")
plot_eff((total_all_seed_hist, bin_edges), (total_matched_seed_hist, bin_edges), "caseed_pur.npz")
