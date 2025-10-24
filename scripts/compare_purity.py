import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("../../fm4npp_eval/script/purity_fm_data_seed1.csv")
pt = data["pT_center"]
purity = data["purity"]
purity_err = data["purity_error"]


plt.figure(figsize=(10, 6))
plt.errorbar(pt,purity,yerr=purity_err,fmt="s", capsize=2, markersize=8,label="Foundation Model", color='red')

plt.legend(fontsize=18)
plt.xlabel('p$_{T}$ [GeV/c]', fontsize=24,labelpad=2)
plt.ylabel('Purity', fontsize=24)
plt.tick_params(axis='both', which='major', labelsize=17)
#plt.title('Matching Efficiency vs. pT')
plt.ylim(0, 1.1)
plt.tick_params(which='both', length=12)
plt.tick_params(direction="in")
plt.tick_params(top=True, right=True)
plt.xlim(0,3)


data_ca = np.load("caseed_pur.npz")
bin_centers_ca = data_ca["bin_centers"]
purity_ca = data_ca["efficiency"]
lower_bounds_ca = data_ca["lower_bounds"]
upper_bounds_ca = data_ca["upper_bounds"]

plt.errorbar(bin_centers_ca, purity_ca, 
             yerr=[lower_bounds_ca, upper_bounds_ca], fmt='o', capsize=2, markersize=8,label='Cellular Automaton',color='#404040')

plt.legend(fontsize=18)
plt.show()