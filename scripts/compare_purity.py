import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

data = pd.read_csv("purity_fm_data_seed1.csv")
pt = data["pT_center"]
purity = data["purity"]
purity_err = data["purity_error"]


df = pd.read_csv("efficiency_fm_data_seed9.csv")

# Extract relevant columns
pte = df["pT_center"]
effe = df["efficiency"]
eff_erre = df["efficiency_error"]

# Plot with error bars



data_ca = np.load("caseed_pur.npz")
bin_centers_ca = data_ca["bin_centers"]
purity_ca = data_ca["efficiency"]
lower_bounds_ca = data_ca["lower_bounds"]
upper_bounds_ca = data_ca["upper_bounds"]
plt.figure(figsize=(10, 6))

plt.errorbar(bin_centers_ca, purity_ca, 
             yerr=[lower_bounds_ca, upper_bounds_ca], fmt='o', capsize=0, markersize=8,label='CA purity', color='#A0A0A0')


dataeff = np.load("caseed_eff.npz")
bin_centers = dataeff["bin_centers"]
efficiency = dataeff["efficiency"]
lower_bounds = dataeff["lower_bounds"]
upper_bounds = dataeff["upper_bounds"]

color = mcolors.to_rgba('gray', alpha=0.5)
plt.errorbar(bin_centers, efficiency, 
            yerr=[lower_bounds, upper_bounds], fmt='s', capsize=0, markersize=8,label='CA efficiency',color='#C0C0C0')
plt.xlabel('p$_{T}$ [GeV/c]', fontsize=24,labelpad=2)
#plt.ylabel('Efficiency', fontsize=24)
plt.tick_params(axis='both', which='major', labelsize=17)
#plt.title('Matching Efficiency vs. pT')
plt.ylim(0, 1.1)
plt.tick_params(which='both', length=12)
plt.tick_params(direction="in")
plt.tick_params(top=True, right=True)
plt.xlim(0,3)
#plt.text(2,0.7, "Cellular Automaton",fontsize=16)
#plt.grid(True, alpha=0.3)


#plt.figure(figsize=(10, 6))
plt.errorbar(pt,purity,yerr=purity_err,fmt="o", capsize=2, markersize=8,label="FM purity",  color='blue')

plt.legend(fontsize=18)
plt.xlabel('p$_{T}$ [GeV/c]', fontsize=24,labelpad=2)
#plt.ylabel('Purity', fontsize=24)
plt.tick_params(axis='both', which='major', labelsize=17)
#plt.title('Matching Efficiency vs. pT')
plt.ylim(0, 1.1)
plt.tick_params(which='both', length=12)
plt.tick_params(direction="in")
plt.tick_params(top=True, right=True)
plt.xlim(0,3)



plt.errorbar(pte, effe, yerr=eff_erre, fmt="s", capsize=2, markersize=8,
            label="FM efficiency",color='red')
    
    

handles, labels = plt.gca().get_legend_handles_labels()

# Reorder: FM purity (idx 2), FM efficiency (idx 3), CA purity (idx 0), CA efficiency (idx 1)
order = [3, 2, 1, 0]
plt.legend([handles[i] for i in order], [labels[i] for i in order], fontsize=18)

plt.show()


