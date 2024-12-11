import numpy as np
from pathlib import Path
import pickle
from scipy.stats import chi2, beta
from scipy.interpolate import CubicSpline#PchipInterpolator
from emsel.emsel_util import bh_correct, get_1d_s_data_from_type
from emsel.core import HMM
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import seaborn as sns
from cycler import cycler

###### MODIFY

data_dir = "data/ibdne"
output_dir = "output/ibdne"

###### DO NOT MODIFY

plt.rcParams.update({'font.size': 9,
                     'text.usetex': False,
                     'font.family': 'serif',
                     'font.serif': 'cmr10',
                     'mathtext.fontset': 'cm',
                     'axes.unicode_minus': False,
                     'axes.formatter.use_mathtext': True,})

plt.rcParams["axes.prop_cycle"] = cycler(color=sns.color_palette("dark"))

hmm_Nes = np.geomspace(500,500000, 21, dtype=int)
prefixes = ["g125_dal_special_realmatch", "g125_dal_special_ibdne"]
EM_dirs = ["EM/real_matched/boxplots", "EM/ibdne/boxplots"]
subfigs = ["A", "B"]
num_batches = 25

spacing_frac = 2

Nes_space = np.geomspace(2000, 100000, 5000)
true_maxes = [9715, 26205] #ibdne estimated as harmonic mean of ibdne_raw.txt

assert len(true_maxes) == len(prefixes)

all_ests = []
for i in range(len(prefixes)):
    all_ests.append(np.zeros(num_batches))


fig, axs = plt.subplots(1,len(true_maxes),figsize=(6.2,3.1),layout="constrained")
fig.get_layout_engine().set(wspace=.08)
for p_i, prefix in enumerate(prefixes):
    #axs2[p_i%2, p_i//2].set_title(prefixes[p_i])
    axs[p_i].set_xscale("log")
    axs[p_i].set_xlabel(r"$N_e$")
    axs[p_i].set_ylabel("Log-likelihood (a.u.)")
    for i in tqdm(range(num_batches)):
        Ne_lls = np.zeros_like(hmm_Nes, dtype=float)
        condo_lls = np.zeros_like(hmm_Nes, dtype=float)
        for Ne_i, hmm_Ne in enumerate(hmm_Nes):
            em_path = Path(f"{EM_dirs[p_i]}/neutral_{prefix}_seed{100+i}_HMMNe{hmm_Ne}_cond_unif_EM.pkl")
            with open(em_path, "rb") as file:
                hf = pickle.load(file)
            uncondo_sum = hf["neutral_ll"].sum()
            condo_sum = uncondo_sum - hf["cond_correction_ll"].sum()
            Ne_lls[Ne_i] += uncondo_sum
            condo_lls[Ne_i] += condo_sum
        unif_spline = CubicSpline(hmm_Nes, Ne_lls)
        unif_spline_output = unif_spline(Nes_space)
        condo_spline = CubicSpline(hmm_Nes, condo_lls)
        condo_spline_output = condo_spline(Nes_space)

        condo_spline_normalized = condo_spline_output-np.max(condo_spline_output)
        if i==0:
            offset_starter = np.min(condo_spline_normalized)
        offset = offset_starter*i*spacing_frac
        temp_line, = axs[p_i].plot(Nes_space, condo_spline_normalized-offset)
        axs[p_i].plot(Nes_space[np.argmax(condo_spline_normalized)], np.max(condo_spline_normalized)-offset, "*", c=temp_line.get_color(), ms=10)

for i in range(len(true_maxes)):
    axs[i].set_yticks([])
    axs[i].text(-.08, .97, rf"$\bf{{{subfigs[i]}}}$", fontsize=13, transform=axs[i].transAxes)

fig.savefig(Path(output_dir)/f"condo_batch_profile_plots.pdf", format="pdf", bbox_inches="tight")
plt.close(fig)