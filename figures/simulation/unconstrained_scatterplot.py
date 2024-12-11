import numpy as np
import pickle
import matplotlib.pyplot as plt
from itertools import product as itprod
from pathlib import Path
import seaborn as sns
import pandas as pd
from cycler import cycler
from emsel.emsel_util import params_dict_to_str, get_1d_s_data_from_type, convert_from_abbrevs

####### MODIFY

sel_strs = [.05]
num_gens_list = [251]
init_dists = [.25]

data_dir = "data/pure_sim"
EM_dir = "EM/pure_sim"
output_dir = "output/pure_sim"
classified_dir = "classified"



####### DO NOT MODIFY

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
colorlist = ["#1D6996", *[coolormap(i) for i in [1,0]], colors[3], colors[4]]
init_colorlist = colorlist
plt.rcParams["axes.prop_cycle"] = cycler(color=colorlist)

sel_types = ["neutral", "add", "dom", "rec", "over", "under"]
cond_types = ["neutral", "add", "dom", "rec"]

sel_coords_1 = [0, 0.025, 0.05, 0, 0.05, -0.05]
sel_coords_2 = [0, 0.05, 0.05, 0.05, 0, 0]

run_types = ["add", "dom", "rec", "het"]
colorlist = init_colorlist
plt.rcParams["axes.prop_cycle"] = cycler(color=init_colorlist)
for n_i, num_gens in enumerate(num_gens_list):
     for d_i, init_dist in enumerate(init_dists):

        fig, axs = plt.subplots(1, 1, figsize=(4.5, 4.5), layout="constrained")
        if init_dist == "recip":
            ic_str = "1/x"
        elif init_dist == "real_special":
            ic_str = "real"
        else:
            ic_str = f"delta (p={init_dist})"
        neutral_list = []
        s_vals = []
        s_types = []
        s_strs = []
        init_conds = []
        freqs_list = []
        min_quantile = 0
        max_quantile = 0
        runtime_total = 0
        i = -1
        for sel_type, sel_str in itprod(sel_types, sel_strs):
            pdict = {"sel_type": sel_type, "num_gens": num_gens, "sel_str": sel_str, "init_dist": init_dist}
            exp_name = params_dict_to_str(**pdict)
            hmm_filename = Path(f"{EM_dir}/{exp_name}_EM.pkl")
            if not hmm_filename.is_file():
                print(f"hmm file not found: {hmm_filename}")
                continue

            with open(hmm_filename, "rb") as file:
                hf = pickle.load(file)

            axs.scatter(hf["full_run"]["s_final"][0,:], hf["full_run"]["s_final"][1,:], color=colorlist[i] if i >= 0 else "black", s=5, label=convert_from_abbrevs(sel_type))
            axs.plot(sel_coords_1[i+1], sel_coords_2[i+1], ls=None, marker="*", ms=10, color=colorlist[i] if i >= 0 else "black", zorder=3, markeredgewidth=.75, markeredgecolor="k" if i>=0 else "r")
            i += 1
            #pf = np.loadtxt(pd_filename, delimiter="\t")
        axs.axhline(lw=1.5, color="k")
        axs.axvline(lw=1.5, color="k")
        axs.set_xlim([-.1,.1])
        axs.set_ylim([-.1,.1])
        axs.set_xlabel(r"$s_1$")
        axs.set_ylabel(r"$s_2$")
        axs.legend(labelspacing=.2, handlelength=1.5, handleheight=.5, handletextpad=.4, borderpad=.2, borderaxespad=.2,)
        plt.savefig(f"{output_dir}/g{num_gens}_d{init_dist}_unconstrained_scatterplot.pdf", format="pdf", bbox_inches="tight")
        plt.close(fig)