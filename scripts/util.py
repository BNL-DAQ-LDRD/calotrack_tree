# %%
import uproot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
from scipy.stats import binom
from scipy.stats import beta
from scipy.optimize import least_squares

def get_clusters(data, event_id):
    cid = data["reco_cluster_id"].array(library="np")[event_id]
    x = data["reco_cluster_x"].array(library="np")[event_id]
    y = data["reco_cluster_y"].array(library="np")[event_id]
    z = data["reco_cluster_z"].array(library="np")[event_id]
    detid = data["reco_cluster_detid"].array(library="np")[event_id]
    ptid = data["reco_cluster_g4hit_trkid"].array(library="np")[event_id]
    df = pd.DataFrame({"cid": cid, "x": x, "y": y, "z": z, "detid": detid, "ptid": ptid})
    return df

# %%
def get_seeds(data, event_id, clusters, cid_to_index):
    sid = data["tpc_seeds_id"].array(library="np")[event_id]
    cids = data["tpc_seeds_clusters"].array(library="np")[event_id]
    ncid = data["tpc_seeds_nclusters"].array(library="np")[event_id]
    cid_start = data["tpc_seeds_start_idx"].array(library="np")[event_id]
    list_of_cids = [cids[start : start + length] 
                for start, length in zip(cid_start, ncid)]
    lists_of_x = [[clusters.iloc[cid_to_index[cid]]['x'] for cid in cids] for cids in list_of_cids]
    lists_of_y = [[clusters.iloc[cid_to_index[cid]]['y'] for cid in cids] for cids in list_of_cids]
    lists_of_z = [[clusters.iloc[cid_to_index[cid]]['z'] for cid in cids] for cids in list_of_cids]
    df = pd.DataFrame({"sid": sid, "cids": list_of_cids, "x": lists_of_x, "y": lists_of_y, "z": lists_of_z})
    return df

# %%
def get_particles(data, event_id, clusters, cid_to_index, pid_to_cids):
    ppid = data["particle_pid"].array(library="np")[event_id]
    ptid = data["particle_track_id"].array(library="np")[event_id]
    cids = [pid_to_cids.get(p, [])  for p in ptid]
    cids = [[cid for cid in particle_cids if cid in cid_to_index and clusters.iloc[cid_to_index[cid]]['detid'] == 2] for particle_cids in cids]
    lists_of_x = [[clusters.iloc[cid_to_index[cid]]['x'] for cid in cids] for cids in cids]
    lists_of_y = [[clusters.iloc[cid_to_index[cid]]['y'] for cid in cids] for cids in cids]
    lists_of_z = [[clusters.iloc[cid_to_index[cid]]['z'] for cid in cids] for cids in cids]
    px = data["particle_px"].array(library="np")[event_id]
    py = data["particle_py"].array(library="np")[event_id]
    pz = data["particle_pz"].array(library="np")[event_id]
    pt = [np.sqrt(px**2 + py**2) for px, py in zip(px, py)]
    eta = np.arctanh(pz / np.sqrt(px**2 + py**2 + pz**2))
    vz = data["particle_vtx_z"].array(library="np")[event_id]
    df = pd.DataFrame({"ptid": ptid, "ppid": ppid, "pt": pt, "pz": pz, "vz": vz, "eta": eta, "cids": cids, "x": lists_of_x, "y": lists_of_y, "z": lists_of_z})
    # df = df[df['ptid'] > 0]
    # df = df[df['pt'] > 0.1]
    # charge_stable_ppids = [13, -13, 11, -11, 211, -211, 321, -321, 2212, -2212]
    # df = df[df['ppid'].isin(charge_stable_ppids)]
    df = df[df['cids'].apply(len) > 2]
    return df

# Helper function to short-circuit intersection counting
def has_ncommon_or_more(set_a, set_b, n):
    """
    Return True if sets `set_a` and `set_b` share at least `n` elements,
    checking incrementally so we can stop early.
    """
    # Always iterate over the smaller set to reduce lookups
    if len(set_a) <= len(set_b):
        smaller, larger = set_a, set_b
    else:
        smaller, larger = set_b, set_a

    count = 0
    for item in smaller:
        if item in larger:
            count += 1
            if count >= n:
                return True
    return False

def match_seeds_to_particles_optimized(seeds, particles, ncommon):
    particle_cid_sets = []
    for _, particle_row in particles.iterrows():
        particle_cid_sets.append( (particle_row['pt'], set(particle_row['cids'])) )
    
    matched_pts = []
    matches = pd.Series(False, index=seeds.index)
    matchedpts = pd.Series(0.0, index=seeds.index)
    matchedtruthlengths = pd.Series(0, index=seeds.index)
    for irow, seed_row in seeds.iterrows():
        seed_cids = set(seed_row['cids'])
        for ptid, particle_set in particle_cid_sets:
            nmajor = 0.5*max(len(seed_cids), len(particle_set))
            
            if has_ncommon_or_more(seed_cids, particle_set, nmajor):
                
                matched_pts.append(seed_row['sid'])
                matches[irow] = True
                matchedpts[irow] = ptid
                matchedtruthlengths[irow] = len(particle_row['cids'])
                break
    seeds['matched'] = matches
    seeds['truthpt'] = matchedpts
    seeds['truthlength'] = matchedtruthlengths
    return matched_pts
# %%
def match_particles_to_seeds_optimized(particles, seeds, ncommon):
    """
    For each particle in the particles DataFrame, find if there is at least one
    seed in the seeds DataFrame sharing at least `ncommon` cids. If found,
    store the particle's pt in a list.

    Parameters:
        particles (pd.DataFrame): Must have columns 'pt' and 'cids'.
        seeds (pd.DataFrame): Must have columns 'sid' and 'cids'.
        ncommon (int): Threshold for minimum shared cids.

    Returns:
        list: List of particle pt's that matched at least one seed.
    """

    # Precompute sets of cids for each seed to avoid repeated set() calls
    # We'll just store a list of (sid, set_of_cids) so we can iterate quickly.
    seed_cid_sets = []
    for _, seed_row in seeds.iterrows():
        seed_cid_sets.append( (seed_row['sid'], set(seed_row['cids'])) )

    matched_pts = []
    matches = pd.Series(False, index=particles.index)



    # For each particle, check if at least one seed matches
    for irow, particle_row in particles.iterrows():
        particle_cids = set(particle_row['cids'])
        for sid, seed_set in seed_cid_sets:
            nmajor = 0.5*max(len(particle_cids), len(seed_set))
            if has_ncommon_or_more(particle_cids, seed_set, nmajor):
                matched_pts.append(particle_row['pt'])
                # print(f"irow {irow}, len(particles) {len(particles)}")
                matches[irow] = True
                break  # no need to check further seeds once matched
    
    # Add the match status as a column to the particles DataFrame
    particles['matched'] = matches

    return matched_pts


# %%

# Suppose you already have these:
# all_pt = particles['pt'].tolist()
# matched_pt = match_particles_to_seeds_optimized(particles, seeds, ncommon)

# def plot_eff(all_pt, matched_pt):
#     # 1) Plot histograms of all_pt and matched_pt
#     bin_edges = np.linspace(0, 2, 21)  # Define bin edges
#     print(f"bin_edges: {bin_edges}")
#     plt.hist(all_pt, bins=bin_edges, alpha=0.5, label='All')
#     plt.hist(matched_pt, bins=bin_edges, alpha=0.5, label='Matched')
#     plt.xlabel('pt')
#     plt.ylabel('Frequency')
#     plt.title('Histogram of All vs. Matched pT')
#     plt.legend()
#     plt.show()

#     # 2) Compute Efficiency per bin
#     bin_counts_all, bin_edges = np.histogram(all_pt, bins=bin_edges)
#     bin_counts_matched, _     = np.histogram(matched_pt, bins=bin_edges)
#     print(f"bin_counts_all: {bin_counts_all}, bin_counts_matched: {bin_counts_matched}")
#     print(f"bin_edges: {bin_edges}")

#     # Avoid divide-by-zero by checking bin_counts_all before dividing
#     efficiency = np.zeros_like(bin_counts_all, dtype=float)
#     mask = bin_counts_all > 0
#     efficiency[mask] = bin_counts_matched[mask] / bin_counts_all[mask]
#     print(f"efficiency: {efficiency}")

#     # 3) Plot Efficiency vs. pT
#     # Use the bin centers for plotting
#     bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

#     plt.plot(bin_centers, efficiency, marker='o', linestyle='-')
#     plt.xlabel('pt')
#     plt.ylabel('Efficiency')
#     plt.title('Matching Efficiency vs. pT')
#     plt.ylim(0, 1)  # Efficiency ranges from 0 to 1
#     plt.show()


# Define circle fit function
def fit_circle(x, y):
    """
    Fit a circle to 2D points using least squares
    Returns (center_x, center_y, radius)
    """
    x = np.array(x)
    y = np.array(y)
    
    def calc_residuals(params, x, y):
        xc, yc, r = params
        return np.sqrt((x - xc)**2 + (y - yc)**2) - r
    
    # Initial guess: center at mean, radius as std
    x_m = np.mean(x)
    y_m = np.mean(y)
    r_guess = np.sqrt((x - x_m)**2 + (y - y_m)**2).mean()
    
    result = least_squares(calc_residuals, [x_m, y_m, r_guess], args=(x, y))
    xc, yc, r = result.x
    return xc, yc, r

# Apply to each row
def fit_circle_row(row):
    x = np.array(row['x'])
    y = np.array(row['y'])
    
    # Check for minimum points
    if len(x) < 3:
        return pd.Series({
            'center_x': np.nan,
            'center_y': np.nan,
            'radius': np.nan,
            'n_points': len(x),
            'fit_status': 'too_few_points'
        })
    
    # Check for NaN
    if np.any(np.isnan(x)) or np.any(np.isnan(y)):
        return pd.Series({
            'center_x': np.nan,
            'center_y': np.nan,
            'radius': np.nan,
            'n_points': len(x),
            'fit_status': 'contains_nan'
        })
    
    try:
        center_x, center_y, radius = fit_circle(x, y)
        return pd.Series({
            'center_x': center_x,
            'center_y': center_y,
            'radius': radius,
            'n_points': len(x),
            'fit_status': 'success'
        })
    except Exception as e:
        return pd.Series({
            'center_x': np.nan,
            'center_y': np.nan,
            'radius': np.nan,
            'n_points': len(x),
            'fit_status': f'fit_failed: {str(e)}'
        })
def plot_eff(all_pt_hist, matched_pt_hist, filename):
    """
    Plot efficiency from histogram data.
    
    Parameters:
        all_pt_hist: Tuple of (counts, bin_edges) from np.histogram
        matched_pt_hist: Tuple of (counts, bin_edges) from np.histogram
    """
    # Extract histogram data
    bin_counts_all, bin_edges = all_pt_hist
    bin_counts_matched, _ = matched_pt_hist
    
    # 1) Plot histograms using the bin data
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # visualize the input histograms
    # plt.figure(figsize=(10, 6))
    # plt.bar(bin_centers, bin_counts_all, width=bin_edges[1]-bin_edges[0], 
    #         alpha=0.5, label='All', align='center')
    # plt.bar(bin_centers, bin_counts_matched, width=bin_edges[1]-bin_edges[0], 
    #         alpha=0.5, label='Matched', align='center')
    # plt.xlabel('pT (GeV/c)')
    # plt.ylabel('Frequency')
    # plt.title('Histogram of All vs. Matched pT')
    # plt.legend()
    # plt.grid(True, alpha=0.3)
    # plt.show()

    # 2) Compute and plot efficiency
    efficiency = np.zeros_like(bin_counts_all, dtype=float)
    mask = bin_counts_all > 0
    efficiency[mask] = bin_counts_matched[mask] / bin_counts_all[mask]
    
    # 3) Compute asymmetric confidence intervals
    confidence_level = 0.68  # ~1σ
    lower_bounds = np.zeros_like(efficiency)
    upper_bounds = np.zeros_like(efficiency)

    # binomial confidence intervals
    for i in range(len(bin_counts_all)):
        n = bin_counts_all[i]
        k = bin_counts_matched[i]
        if n > 0:
            ci_low, ci_up = binom.interval(confidence_level, n, k / n)
            lower_bounds[i] = efficiency[i] - ci_low / n
            upper_bounds[i] = ci_up / n - efficiency[i]
        else:
            lower_bounds[i] = 0
            upper_bounds[i] = 0

    # beta distribution confidence intervals
    # for i in range(len(bin_counts_all)):
    #     n = bin_counts_all[i]
    #     k = bin_counts_matched[i]
    #     if n > 0:
    #         alpha_post = k + 0.5
    #         beta_post = n - k + 0.5
    #         ci_low = beta.ppf((1 - confidence_level) / 2, alpha_post, beta_post)
    #         ci_up  = beta.ppf(1 - (1 - confidence_level) / 2, alpha_post, beta_post)
    #         lower_bounds[i] = efficiency[i] - ci_low
    #         upper_bounds[i] = ci_up - efficiency[i]

    # 4) Plot efficiency with asymmetric error bars
    plt.figure(figsize=(10, 6))
    plt.errorbar(bin_centers, efficiency, 
                 yerr=[lower_bounds, upper_bounds], fmt='o', capsize=2, markersize=7)
    plt.xlabel('p$_{T}$ [GeV/c]', fontsize=17)
    plt.ylabel('Efficiency', fontsize=17)
    plt.tick_params(axis='both', which='major', labelsize=14)
    #plt.title('Matching Efficiency vs. pT')
    plt.ylim(0, 1.1)
    plt.tick_params(which='both', length=10)
    plt.tick_params(direction="in")
    plt.tick_params(top=True, right=True)
    plt.xlim(0,3)
    #plt.grid(True, alpha=0.3)
    plt.show()

    np.savez(filename,
         bin_centers=bin_centers,
         efficiency=efficiency,
         lower_bounds=lower_bounds,
         upper_bounds=upper_bounds)
    # Print some statistics
    print(f"Overall efficiency: {sum(bin_counts_matched)/sum(bin_counts_all):.4f}")
    return efficiency, bin_centers

def get_group_ids(groups, n, cid_to_index, cid_to_pt, bin_edges):
    """
    Maps each point to its group ID.
    
    Parameters:
        groups is a pandas DataFrame with a column 'cids' containing lists of point indices
        n (int): Total number of points
        cid_to_pt (dict): Mapping from cid to pt
        bin_edges (list): Bin edges for histogramming
        
    Returns:
        list of list, first index is corresponding to bin index
    """
    # Initialize a list of lists to store the ids for each pt bin
    bin_ids = [[-1]*n for _ in range(len(bin_edges) - 1)]
    
    # Process each group
    for groupid, group in groups.iterrows():
        for cid in group['cids']:
            if cid in cid_to_pt:
                pt = cid_to_pt[cid]
                # Find which bin this pt belongs to
                bin_idx = np.digitize(pt, bin_edges) - 1
                # Only store if it falls within our bins
                if 0 <= bin_idx < len(bin_edges) - 1:
                    # If cid has an index, add it to the appropriate bin
                    if cid in cid_to_index:
                        bin_ids[bin_idx][cid_to_index[cid]] = groupid # overwrite if already exists
    return bin_ids

def cid_to_pt_mapping(particles):
    """
    Create a mapping from cid to pt.
    
    Parameters:
        particles is a pandas DataFrame with a column 'cids' containing lists of point indices
        n (int): Total number of points
        
    Returns:
        dict: Mapping from cid to pt
    """
    cid_to_pt = {}
    for _, particle in particles.iterrows():
        for cid in particle['cids']:
            cid_to_pt[cid] = particle['pt'] # overwrite if already exists
    return cid_to_pt
