import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

cadata = np.load("ntpc_caseed_data.npz")
fmdata = np.load("../../fm4npp_eval/script/ntpc_fm_seed1.npz")
fmcounts = fmdata['counts']
fmbin_edges = fmdata['bin_edges']
normalized_fmcounts = fmcounts / fmcounts.sum()
normalized_fmcounts = fmcounts

cadata = np.load('ntpc_caseed_data.npz')
cacounts = cadata['counts']
normalized_cacounts = cacounts / cacounts.sum()
normalized_cacounts = cacounts
cabin_edges = cadata['bin_edges']
print("ca counts sum " + str(cacounts.sum()))
print("fm counts sum " + str(fmcounts.sum()))

plt.bar(cabin_edges[:-1], cacounts, width=np.diff(cabin_edges), align='edge', edgecolor='black')
plt.xlabel('nTPC')
plt.ylabel('Count')
plt.show()

# Method 2: Plot as step histogram (looks more like original hist)
bin_centers = (cabin_edges[:-1] + cabin_edges[1:]) / 2
plt.step(cabin_edges[:-1], normalized_cacounts, where='post',label='Cellular Automaton', color='#404040')
plt.xlabel('nTPC')
plt.xlim(20,50)
plt.ylabel('Arbitrary')

#plt.show()

fmbincenters = (fmbin_edges[:-1] + fmbin_edges[1:]) / 2
plt.step(fmbin_edges[:-1], normalized_fmcounts, where='post',label="Foundation Model",color='red')
plt.legend(fontsize=15)
plt.show()
