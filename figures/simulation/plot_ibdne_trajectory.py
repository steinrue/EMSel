import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pickle
from scipy.stats import chi2

data_dir = "data"
output_dir = "output/ibdne"
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

def cmap(val):
    return (1-val, 0, val, 1)
coolormap = plt.get_cmap("Dark2")
colorlist = ["#1D6996", *[coolormap(i) for i in [1,0]], colors[3], colors[4], colors[5]]


fig, axs = plt.subplots(1,1,figsize=(3.6,3.1), layout="constrained")
ibdne_traj = np.loadtxt(f"{data_dir}/ibdne_raw.txt")
ibdne_gens = np.loadtxt(f"{data_dir}/ibdne_original.txt")[:, 0][::-1]
axs.plot(ibdne_gens, ibdne_traj, lw=1.75)
axs.set_xlabel("Generations before present")
axs.set_ylabel(r"$N_e$")
axs.set_xlim([0, 200])
axs.set_yscale("log")
axs.set_ylim([1e3, 1e8])
axs.grid(True)
axs.axvline(ibdne_gens[0], color="k", ls="--", lw=1)
axs.axvline(ibdne_gens[-1], color="k", ls="--", lw=1)
fig.savefig(Path(f"{output_dir}/ibdne_traj.pdf"), format="pdf", bbox_inches="tight")
plt.close(fig)