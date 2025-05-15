import numpy as np
from pathlib import Path
import pickle
from scipy.stats import chi2, beta
from scipy.interpolate import CubicSpline
from emsel.emsel_util import bh_correct, get_1d_s_data_from_type
from emsel.core import HMM
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.ticker import ScalarFormatter
from copy import deepcopy

###### MODIFY

data_dir = "data/pure_sim/boxplots"
EM_dir = "EM/pure_sim/boxplots"
output_dir = "output/pure_sim"

# sim_Nes = [2500, 10000, 40000]
# prefixes = [*[f"g251_d25_Ne{sim_Ne}" for sim_Ne in sim_Nes]]
# labels = [rf"$N_e = {{{str(i)}}}$" for i in sim_Nes]
# EM_dirs = ["EM/pure_sim/boxplots", "EM/pure_sim/boxplots", "EM/pure_sim/boxplots"]
# true_maxes = [2500, 10000, 40000]

prefixes = ["g125_dal_special_realmatch", "g125_dal_special_ibdne"]#
labels = ["Real matched", "IBDNe"]#
EM_dirs = ["EM/real_matched/boxplots", "EM/ibdne/boxplots"]#
true_maxes = [9715,26205] #ibdne estimated as harmonic mean of ibdne_raw.txt

###### DO NOT MODIFY

plt.rcParams.update({'font.size': 9,
                     'text.usetex': False,
                     'font.family': 'serif',
                     'font.serif': 'cmr10',
                     'mathtext.fontset': 'cm',
                     'axes.unicode_minus': False,
                     'axes.formatter.use_mathtext': True,})

hmm_Nes = np.geomspace(500,500000, 21, dtype=int)


run_type = "unif"
num_batches = 25

Nes_space = np.geomspace(2000, 100000, 5000)




assert len(true_maxes) == len(prefixes)
tick_labels = deepcopy(true_maxes)

all_ests = []
all_labels = []
for p_i, prefix in enumerate(prefixes):
    for i in tqdm(range(num_batches)):
        Ne_lls = np.zeros_like(hmm_Nes, dtype=float)
        condo_lls = np.zeros_like(hmm_Nes, dtype=float)
        for Ne_i, hmm_Ne in enumerate(hmm_Nes):
            seed = 100+i
            em_path = Path(f"{EM_dirs[p_i]}/neutral_{prefix}_seed{seed}_HMMNe{hmm_Ne}_cond_{run_type}_EM.pkl")
            with open(em_path, "rb") as file:
                hf = pickle.load(file)
            uncondo_sum = hf["neutral_ll"].sum()
            condo_sum = uncondo_sum - hf["cond_correction_ll"].sum()
            Ne_lls[Ne_i] += uncondo_sum
            condo_lls[Ne_i] += condo_sum


            # if Ne_i == 0 and chrom == 1:
            #     print(hf["neutral_itercount"])
            # Ne_lls[Ne_i] += hf["neutral_ll"].sum()
        unif_spline = CubicSpline(hmm_Nes, Ne_lls)
        unif_spline_output = unif_spline(Nes_space)
        condo_spline = CubicSpline(hmm_Nes, condo_lls)
        condo_spline_output = condo_spline(Nes_space)

        all_ests.append(Nes_space[np.argmax(condo_spline_output)])
        all_labels.append(labels[p_i])

#print([f"{est:.0f}" for est in all_ests])

massaged_data = zip(all_ests, all_labels)
data_df = pd.DataFrame(massaged_data, columns=[r"$\hat{N}_e$", "Dataset"])

box_lws = .7
fig, axs = plt.subplots(1,1,figsize=(3.1,3.1),layout="constrained")
axs.set_yscale("log")
axs.set_yticks(true_maxes)
axs.set_yticklabels(tick_labels)
for true_max in true_maxes:
    axs.axhline(true_max, color="r", ls="--", lw=1)
sns.boxplot(data=data_df, y=r"$\hat{N}_e$", x="Dataset", hue="Dataset", ax=axs, widths=.4, boxprops={"lw": box_lws},
            medianprops={"lw": box_lws}, whiskerprops={"lw": box_lws},
            capprops={"lw": box_lws}, flierprops={"alpha": .7, "ms": 5})
axs.set_yticks(true_maxes)
axs.get_yaxis().set_major_formatter(ScalarFormatter())
axs.set_ylabel(r"$\hat{N}_e$")
axs.set_xlabel("Dataset")

if "IBDNe" in labels:
    ibdneloc = labels.index("IBDNe")
    axs.get_yticklabels()[ibdneloc].set_color("red")
fig.savefig(Path(output_dir)/f"condo_batch_boxplots_final.pdf", format="pdf", bbox_inches="tight")
plt.close(fig)