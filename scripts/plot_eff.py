import pandas as pd
from util import *

# Open the HDF5 file in read mode
with pd.HDFStore('data_event_0.h5', mode='r') as store:
    clusters = store['clusters']      # Saved as table format
    seeds = store['seeds']            # Saved as fixed format
    particles = store['particles']    # Saved as fixed format
    

    ncommon = 5  # change to your desired threshold
    matched_pt = match_particles_to_seeds_optimized(particles, seeds, ncommon)
    all_pt = particles['pt'].tolist()
    print(f"All {len(all_pt)} matched: {len(matched_pt)}")

    plot_eff(all_pt, matched_pt)