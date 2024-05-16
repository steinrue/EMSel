import numpy as np
import pickle
import matplotlib.pyplot as plt
from itertools import product as itprod
from pathlib import Path
import seaborn as sns
import pandas as pd
from cycler import cycler
from emsel.emsel_util import params_dict_to_str, get_1d_s_data_from_type, convert_from_abbrevs

###### MODIFY

data_dir = "data"
EM_dir = "EM"
output_dir = "output"

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
colorlist = ["#1D6996", *[coolormap(i) for i in [1,0]], colors[3], colors[4]]
init_colorlist = colorlist
plt.rcParams["axes.prop_cycle"] = cycler(color=colorlist)

sel_types = ["neutral", "add", "dom", "rec", "over", "under"]
sel_strs = [.005, .01, .025, .05]
num_gens = 251
init_dist = .25

cond_only = False

fake_sel_types = ["neutral"]
#usable_sel_types = ["neutral", "add", "over", "under"]
usable_types = ["neutral","add", "dom", "rec", "over", "under"]
oned_types = ["add", "dom", "rec", "over", "under"]
run_types = ["add", "dom", "rec", "het"]
cond_types = ["neutral", "add", "dom", "rec"]
iter_types = cond_types if cond_only else sel_types


add_em_fig, add_em_axs = plt.subplots(1, 1, figsize=(3.1, 3.1), layout="constrained")

neutral_list = []
s_vals = []
s_types = []
s_strs = []
init_conds = []
freqs_list = []
min_quantile = 0
max_quantile = 0
runtime_total = 0
illegal_s = 0
for sel_type, sel_str in itprod(iter_types, sel_strs):
    pdict = {"sel_type": sel_type, "num_gens": num_gens, "sel_str": sel_str, "init_dist": init_dist}
    exp_name = params_dict_to_str(**pdict)
    if sel_type == "neutral":
        if exp_name not in neutral_list:
            neutral_list.append(exp_name)
            sel_str = 0.0
        else:
            continue

    hmm_filename = Path(f"{EM_dir}/{exp_name}_EM.pkl")

    pd_filename = Path(f"{data_dir}/{exp_name}_data.csv")
    pdata_filename = Path(f"{data_dir}/{exp_name}_pd.pkl")

    if not pd_filename.is_file():
        print(f"pd file not found: {pd_filename}")
        continue
    if not hmm_filename.is_file():
        print(f"hmm file not found: {hmm_filename}")
        continue


    with open(hmm_filename, "rb") as file:
        hf = pickle.load(file)

    pf = np.loadtxt(pd_filename, delimiter="\t")

    with open(pdata_filename, "rb") as file:
        pdict = pickle.load(file)

    num_pts = hf["neutral_ll"].shape[0]
    idx_list = np.arange(num_pts)

    if sel_type == "neutral":
        s_data = get_1d_s_data_from_type(hf[f"add_run"]["s_final"][:, idx_list[:num_pts]], "add")
        if s_data.shape[0] > 0:
            max_quantile = max(max_quantile, np.quantile(s_data, .99))
            min_quantile = min(min_quantile, np.quantile(s_data, .01))
        s_vals.extend(s_data.tolist())
        s_types.extend([f"add"] * num_pts)
        s_strs.extend(["Neutral"] * num_pts)
    else:
        if sel_type != "under":
            s_data = get_1d_s_data_from_type(hf[f"add_run"]["s_final"][:, idx_list[:num_pts]], "add")
        else:
            s_data = get_1d_s_data_from_type(-hf[f"add_run"]["s_final"][:, idx_list[:num_pts]], "add")
        if s_data.shape[0] > 0:
            max_quantile = max(max_quantile, np.quantile(s_data, .99))
            min_quantile = min(min_quantile, np.quantile(s_data, .01))
        s_vals.extend(s_data.tolist())
        s_types.extend([f"{sel_type}"] * num_pts)
        s_strs.extend([sel_str] * num_pts)
s_types = convert_from_abbrevs(s_types, shortall=True)

massaged_data = zip(s_vals, s_strs, s_types)#

s_stuff = pd.DataFrame(data=massaged_data, columns=[r"$\hat{s}$", r"True $s$", "Mode of selection"])
box = sns.boxplot(data=s_stuff, x=r"True $s$", y=r"$\hat{s}$", hue="Mode of selection", dodge=True, width=.75,
                  ax=add_em_axs,fliersize=1, boxprops={"lw": .5}, medianprops={"lw": .5}, whiskerprops={"lw": .5},
                  capprops={"lw": .5}, flierprops={"alpha": .7})

counts = []
x_coords = []
for s_i, strength in enumerate(["Neutral", 0.005, 0.01, 0.025, 0.05]):
    for m_i, mode in enumerate(["Add.", "Dom.", "Rec.", "Over.", "Under."]):
        bin_mask = ((s_stuff[r"True $s$"]==strength)&(s_stuff["Mode of selection"]==mode)&((s_stuff[r"$\hat{s}$"]>max_quantile)|(s_stuff[r"$\hat{s}$"]<min_quantile)))
        if (bin_counts := bin_mask.sum()) > 0:
            counts.append(bin_counts)
        else:
            counts.append("")
        x_coords.append(s_i+(m_i-2)*.15)
for c_i, count in enumerate(counts):
    add_em_axs.text(x_coords[c_i], max_quantile*1.01, str(count), ha="center", va="bottom", color=colorlist[c_i%len(colorlist)], size=7)

add_em_axs.set_ylim([min_quantile, max_quantile])
add_em_axs.axhline(0, color="r", alpha=.65, ls="--", lw=.5)
for sel_str in sel_strs:
    add_em_axs.axhline(sel_str, color="r", alpha=.6, ls="--", lw=.5)
add_em_axs.text(-.24, 1.01, r"$\bf{A}$", fontsize=13, transform=add_em_axs.transAxes)
add_em_axs.legend(loc="upper center", fontsize=7, labelspacing=.2, handlelength=1.5, handleheight=.5, handletextpad=.4, borderpad=.2, borderaxespad=.2)
plt.savefig(f"{output_dir}/all_add_em_boxplots.pdf", format="pdf", bbox_inches="tight")
plt.close(add_em_fig)

add_sim_fig, add_sim_axs = plt.subplots(1, 1, figsize=(3.1, 3.1), layout="constrained")
neutral_list = []
s_vals = []
s_types = []
s_strs = []
init_conds = []
freqs_list = []
min_quantile = 0
max_quantile = 0
runtime_total = 0
illegal_s = 0
for sel_type, sel_str in itprod(["neutral", "add"], sel_strs):
    pdict = {"sel_type": sel_type, "num_gens": num_gens, "sel_str": sel_str, "init_dist": init_dist}
    exp_name = params_dict_to_str(**pdict)
    if sel_type == "neutral":
        if exp_name not in neutral_list:
            neutral_list.append(exp_name)
            sel_str = 0.0
        else:
            continue

    hmm_filename = Path(f"{EM_dir}/{exp_name}_EM.pkl")

    pd_filename = Path(f"{data_dir}/{exp_name}_data.csv")
    pdata_filename = Path(f"{data_dir}/{exp_name}_pd.pkl")

    if not pd_filename.is_file():
        print(f"pd file not found: {pd_filename}")
        continue
    if not hmm_filename.is_file():
        print(f"hmm file not found: {hmm_filename}")
        continue


    with open(hmm_filename, "rb") as file:
        hf = pickle.load(file)

    pf = np.loadtxt(pd_filename, delimiter="\t")

    with open(pdata_filename, "rb") as file:
        pdict = pickle.load(file)

    num_pts = hf["neutral_ll"].shape[0]
    idx_list = np.arange(num_pts)
    if sel_type == "neutral":
        for run_type in run_types:
            if run_type != "het":
                run_type_val = run_type
            else:
                run_type_val = "over"
            s_data = get_1d_s_data_from_type(hf[f"{run_type}_run"]["s_final"][:, idx_list[:num_pts]], run_type)
            if s_data.shape[0] > 0:
                max_quantile = max(max_quantile, np.quantile(s_data, .99))
                min_quantile = min(min_quantile, np.quantile(s_data, .01))
            s_vals.extend(s_data.tolist())
            s_types.extend([f"{run_type}"] * num_pts)
            s_strs.extend(["Neutral"] * num_pts)
    else:
        for run_type in run_types:
            s_data = get_1d_s_data_from_type(hf[f"{run_type}_run"]["s_final"][:, idx_list[:num_pts]], run_type)
            if s_data.shape[0] > 0:
                max_quantile = max(max_quantile, np.quantile(s_data, .99))
                min_quantile = min(min_quantile, np.quantile(s_data, .01))
            s_vals.extend(s_data.tolist())
            s_types.extend([f"{run_type}"] * num_pts)
            s_strs.extend([sel_str] * num_pts)

s_types = convert_from_abbrevs(s_types, shortall=True)

massaged_data = zip(s_vals, s_strs, s_types)#

s_stuff = pd.DataFrame(data=massaged_data, columns=[r"$\hat{s}$", r"True $s$", "Mode of selection"])
box = sns.boxplot(data=s_stuff, x=r"True $s$", y=r"$\hat{s}$", hue="Mode of selection", dodge=True, width=.75,
                  ax=add_sim_axs,fliersize=1, boxprops={"lw": .5}, medianprops={"lw": .5}, whiskerprops={"lw": .5},
                  capprops={"lw": .5}, flierprops={"alpha": .7})

counts = []
x_coords = []
for s_i, strength in enumerate(["Neutral", 0.005, 0.01, 0.025, 0.05]):
    for m_i, mode in enumerate(["Add.", "Dom.", "Rec.", "Over.", "Under."]):
        bin_mask = ((s_stuff[r"True $s$"]==strength)&(s_stuff["Mode of selection"]==mode)&((s_stuff[r"$\hat{s}$"]>max_quantile)|(s_stuff[r"$\hat{s}$"]<min_quantile)))
        if (bin_counts := bin_mask.sum()) > 0:
            counts.append(bin_counts)
        else:
            counts.append("")
        x_coords.append(s_i+(m_i-2)*.1)
for c_i, count in enumerate(counts):
    add_sim_axs.text(x_coords[c_i], max_quantile*1.01, str(count), ha="center", va="bottom", color=colorlist[c_i%len(colorlist)], size=7)

add_sim_axs.set_ylim([min_quantile, max_quantile])
add_sim_axs.axhline(0, color="r", alpha=.65, ls="--", lw=.5)
for sel_str in sel_strs:
    add_sim_axs.axhline(sel_str, color="r", alpha=.6, ls="--", lw=.5)
add_sim_axs.text(-.24, 1.01, r"$\bf{B}$", fontsize=13, transform=add_sim_axs.transAxes)
add_sim_axs.legend(loc="upper center", fontsize=7, labelspacing=.2, handlelength=1.5, handleheight=.5, handletextpad=.4, borderpad=.2, borderaxespad=.2)
plt.savefig(f"{output_dir}/all_add_sim_boxplots.pdf", format="pdf", bbox_inches="tight")
plt.close(add_sim_fig)