import matplotlib.pyplot as plt
from pathlib import Path
import pickle
from cycler import cycler
import numpy as np

###### MODIFY

output_dir = "output"
genodata_type = "capture_only"

###### DO NOT MODIFY

plt.rcParams.update({'font.size': 9,
                     'text.usetex': False,
                     'font.family': 'serif',
                     'font.serif': 'cmr10',
                     'mathtext.fontset': 'cm',
                     'axes.unicode_minus': False,
                     'axes.formatter.use_mathtext': True,})
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
coolormap = plt.get_cmap("Dark2")
colorlist = ["#1D6996", *[coolormap(i) for i in [1,0]], colors[3], colors[4], colors[5]]
plt.rcParams["axes.prop_cycle"] = cycler(color=colorlist)

with open(Path(f"{output_dir}/GB_v54.1_{genodata_type}_agg_data.pkl"),"rb") as file:
    agg_data = pickle.load(file)

fig, axs = plt.subplots(1,1,figsize=(5,5),layout="constrained")
axs.hist(agg_data["all_missingness"], bins=25)
axs.set_xlim([0,1])
axs.set_xlabel("Missingness fraction")
axs.set_ylabel("Density")
fig.savefig(f"{output_dir}/{genodata_type}_filtered_missingness.pdf", bbox_inches="tight")

fig2, axs2 = plt.subplots(1,1,figsize=(5,5),layout="constrained")
axs2.hist(agg_data["all_means"], bins=np.linspace(0,1,26))
axs2.set_xlim([0,1])
axs2.set_ylim([0,70000])
axs2.set_xlabel("Allele frequency")
axs2.set_ylabel("Density")
fig2.savefig(f"{output_dir}/{genodata_type}_filtered_means_equalspace.pdf", bbox_inches="tight")

# add_means = np.loadtxt(f"{output_dir}/GB_v54.1_capture_only_add_means.txt")
# fig3, axs3 = plt.subplots(1,1,figsize=(5,5),layout="constrained")
# axs3.hist(add_means[1:], bins=np.linspace(0,1,26))
# axs3.set_xlim([0,1])
# axs3.set_ylim([0,70000])
# axs3.set_xlabel("Allele frequency")
# axs3.set_ylabel("Density")
#
# fig3.savefig(f"{output_dir}/{genodata_type}_filtered_add_means_equalspace.pdf", bbox_inches="tight")

