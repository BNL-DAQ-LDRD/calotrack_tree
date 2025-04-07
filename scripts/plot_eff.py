import pandas as pd
import numpy as np
from util import *

bin_edges = np.linspace(0, 2, 21)  # Define bin edges
total_matched_hist = np.zeros(len(bin_edges) - 1, dtype=int)
total_all_hist = np.zeros(len(bin_edges) - 1, dtype=int)

# Loop to read and process files from event_0 to event_10
for ievent in range(1):
    with pd.HDFStore(f'data_event_{ievent}.h5', mode='r') as store:
        clusters = store['clusters']      # Saved as table format
        seeds = store['seeds']            # Saved as fixed format
        particles = store['particles']    # Saved as fixed format

        ncommon = 5  # change to your desired threshold
        matched_pt = match_particles_to_seeds_optimized(particles, seeds, ncommon)
        all_pt = particles['pt'].tolist()

        matched_hist, _ = np.histogram(matched_pt, bins=bin_edges)
        all_hist, _ = np.histogram(all_pt, bins=bin_edges)
        total_matched_hist += matched_hist
        total_all_hist += all_hist
        print(f"Event {ievent}: All {len(all_pt)} matched: {len(matched_pt)}\n")

print("Total Matched histogram:", total_matched_hist)
print("Total All histogram:", total_all_hist)

# Plotting the histograms
plot_eff((total_all_hist, bin_edges), (total_matched_hist, bin_edges))