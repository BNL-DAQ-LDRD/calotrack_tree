import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

data = np.load("caseed_eff.npz")
bin_centers = data["bin_centers"]
efficiency = data["efficiency"]
lower_bounds = data["lower_bounds"]
upper_bounds = data["upper_bounds"]

plt.figure(figsize=(10, 6))
plt.errorbar(bin_centers, efficiency, 
             yerr=[lower_bounds, upper_bounds], fmt='o', capsize=2, markersize=8,label='Cellular Automaton',color='#404040')
plt.xlabel('p$_{T}$ [GeV/c]', fontsize=24,labelpad=2)
plt.ylabel('Efficiency', fontsize=24)
plt.tick_params(axis='both', which='major', labelsize=17)
#plt.title('Matching Efficiency vs. pT')
plt.ylim(0, 1.1)
plt.tick_params(which='both', length=12)
plt.tick_params(direction="in")
plt.tick_params(top=True, right=True)
plt.xlim(0,3)
#plt.text(2,0.7, "Cellular Automaton",fontsize=16)
#plt.grid(True, alpha=0.3)


df = pd.read_csv("efficiency_fm_data_seed9.csv")

# Extract relevant columns
pt = df["pT_center"]
eff = df["efficiency"]
eff_err = df["efficiency_error"]

# Plot with error bars


plt.errorbar(pt, eff, yerr=eff_err, fmt="s", capsize=2, markersize=8,
             label="Foundation Model",color='red')
plt.legend(fontsize=18)

plt.show()

