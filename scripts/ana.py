import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from util import *
from sklearn.metrics import adjusted_rand_score
from util import *


var_name = 'eta'
bin_edges = np.linspace(-2, 2, 40)

var_name = 'vz'
bin_edges = np.linspace(-30, 30, 30)

var_name = 'nclus'
bin_edges = np.linspace(0, 80, 81)


matched_vars =[]
all_vars = []
for ievent in range(0, 100):
    print (f"Processing event {ievent}...")
    # fname = f'2025-05-13-pp-1k-cut/data_event_{ievent}.h5'
    fname = f'2025-05-20-pp-1k-nocut/data_event_{ievent}.h5'
    #fname = f'data_event_{ievent}.h5'
    with pd.HDFStore(fname, mode='r') as store:
        clusters = store['clusters']      # Saved as fixed format
        cid_to_index = {cid: index for index, cid in enumerate(clusters['cid'])}
        seeds = store['seeds']            # Saved as fixed format
        particles = store['particles']    # Saved as fixed format
        particles['nclus'] = particles['cids'].apply(len) # Add nclus column to particles

        ncommon = 30  # change to your desired threshold
        matched_pt = match_particles_to_seeds_optimized(particles, seeds, ncommon)
        matched_particles = particles[particles['matched'] == True]

        matched_vars.extend(particles[particles['matched'] == True][var_name].tolist())
        all_vars.extend(particles[var_name].tolist())
        # print(f"Matched particles: {matched_vars}")
        # print(f"Unmatched particles: {unmatched_vars}")

        doplotting = False
        if doplotting:
            # Create a 2D YX scatter plot for clusters
            plt.figure(figsize=(10, 8))
            # Plot all clusters from particles as transparent gray dots
            all_particle_cids = [cid for p_cids in particles['cids'] for cid in p_cids]
            particle_cluster_indices = [cid_to_index[cid] for cid in all_particle_cids if cid in cid_to_index]
            particle_clusters = clusters.iloc[particle_cluster_indices]
            plt.scatter(particle_clusters['x'], particle_clusters['y'], color='gray', alpha=0.3, s=10, label='Particle Clusters')
            # Plot seeds' clusters with different colors for each seed
            for i, seed in seeds.iterrows():
                seed_cids = seed['cids']
                seed_cluster_indices = [cid_to_index[cid] for cid in seed_cids if cid in cid_to_index]
                if seed_cluster_indices:
                    seed_clusters = clusters.iloc[seed_cluster_indices]
                plt.scatter(seed_clusters['x'], seed_clusters['y'], s=20, label=f'Seed {i}' if i < 10 else None)  # Only show first 10 seeds in legend
            plt.xlabel('X Position')
            plt.ylabel('Y Position')
            plt.title(f'Cluster Distribution in Y-X Plane, Event {ievent}')
            plt.grid(True, alpha=0.3)
            plt.legend(loc='upper right', bbox_to_anchor=(1.15, 1))
            plt.tight_layout()
            plt.savefig(f'cluster_yx_scatter_event_{ievent}.png')
            plt.show()

# Create a single plot for both matched and all variables
plt.figure(figsize=(10, 6))

# Plot histograms on the same axes
plt.hist(all_vars, bins=bin_edges, histtype='step', linewidth=2, color='red', label='All')
plt.hist(matched_vars, bins=bin_edges, histtype='step', linewidth=2, color='blue', label='Matched')

plt.xlabel(f'{var_name}')
plt.ylabel('Counts')
plt.title(f'{var_name} Distribution - Matched vs All')
plt.grid(True, alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig(f'hist_{var_name}.png')
plt.show()



