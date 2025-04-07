# %%
import uproot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
from util import *

# data_file = "../macro/2025-04-01-meeting/testout.root"
# data_file = "../macro/testout-low-centrality-1event.root"
data_file = "../macro/2025-04-06-40event/testout.root"
data = uproot.open(data_file)["T;1"]

interactive_validation = False
if interactive_validation:
    event_id = 0  # pick different events
    with uproot.open(data_file) as fp:
        print(list(fp['T;1'].keys()))

    clusters = get_clusters(data, event_id)
    cid_to_index = {cid: index for index, cid in enumerate(clusters['cid'])}
    pid_to_cids = clusters.groupby('ptid')['cid'].apply(list).to_dict()

    seeds = get_seeds(data, event_id, clusters, cid_to_index)
    print(f"Number of seeds: {len(seeds)}")
    print(seeds.iloc[0])

    print(f"Number of clusters: {len(clusters)}")
    print(list(pid_to_cids.items())[0])
    for i in range(10):
        print(list(pid_to_cids.items())[i])

    clusters[clusters["ptid"] < 1e8]["ptid"].plot(kind="hist")

    particles = get_particles(data, event_id, clusters, cid_to_index, pid_to_cids)

    print(f"Number of particles: {len(particles)}")
    for i in range(len(particles)):
        print(f"ptid: {particles.iloc[i]['ptid']}, pt: {particles.iloc[i]['pt']}, pz: {particles.iloc[i]['pz']}, ppid: {particles.iloc[i]['ppid']}, cids: {len(particles['cids'].iloc[i])}")

    ncommon = 5  # change to your desired threshold
    matched_pt = match_particles_to_seeds_optimized(particles, seeds, ncommon)
    all_pt = particles['pt'].tolist()
    print(f"All {len(all_pt)} matched: {len(matched_pt)}")

# %%

for ievent in range(1):
    clusters = get_clusters(data, ievent)
    cid_to_index = {cid: index for index, cid in enumerate(clusters['cid'])}
    pid_to_cids = clusters.groupby('ptid')['cid'].apply(list).to_dict()

    seeds = get_seeds(data, ievent, clusters, cid_to_index)
    print(f"Number of seeds: {len(seeds)}")
    print(seeds.iloc[0])

    print(f"Number of clusters: {len(clusters)}")
    print(list(pid_to_cids.items())[0])
    for i in range(10):
        print(list(pid_to_cids.items())[i])

    clusters[clusters["ptid"] < 1e8]["ptid"].plot(kind="hist")

    particles = get_particles(data, ievent, clusters, cid_to_index, pid_to_cids)

    print(f"Number of particles: {len(particles)}")
    for i in range(len(particles)):
        print(f"ptid: {particles.iloc[i]['ptid']}, pt: {particles.iloc[i]['pt']}, pz: {particles.iloc[i]['pz']}, ppid: {particles.iloc[i]['ppid']}, cids: {len(particles['cids'].iloc[i])}")

    with pd.HDFStore(f'data_event_{ievent}.h5', mode='w', complevel=9, complib='blosc') as store:
        store.put('clusters', clusters, format='table')   # if all scalar columns
        store.put('seeds', seeds, format='fixed')         # workaround for list columns
        store.put('particles', particles, format='fixed') # workaround for list columns


