import numpy as np
import pickle
import matplotlib.pyplot as plt
from itertools import product as itprod
from pathlib import Path
import seaborn as sns
import pandas as pd
from cycler import cycler
from emsel_util import params_dict_to_str, get_1d_s_data_from_type, convert_from_abbrevs

####### MODIFY

sel_strs = [.005, .01, .025, .05]
num_gens_list = [101, 251, 1001]
init_dists = [.005, .25, "recip"]

#set to True for strip plots (e.g. Figure 6)
cond_only = False

data_dir = "data"
EM_dir = "EM"
output_dir = "output"
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

run_types = ["add", "dom", "rec", "het"]
iter_types = cond_types if cond_only else sel_types

for n_i, num_gens in enumerate(num_gens_list):
     for d_i, init_dist in enumerate(init_dists):
        fig, axs = plt.subplots(1, 1, figsize=(3.1, 3.1), layout="constrained")
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
        illegal_s = 0
        for sel_type, sel_str in itprod(iter_types, sel_strs):
            if init_dist in [.005, "recip"] and sel_type == "rec":
                continue
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
            bh_filename = Path(f"{classified_dir}/{exp_name}_classified.pkl")

            if not pd_filename.is_file():
                print(f"pd file not found: {pd_filename}")
                continue
            if not hmm_filename.is_file():
                print(f"hmm file not found: {hmm_filename}")
                continue
            if (cond_only and not bh_filename.is_file()):
                print(f"bh file not found: {bh_filename}")
                continue


            with open(hmm_filename, "rb") as file:
                hf = pickle.load(file)

            pf = np.loadtxt(pd_filename, delimiter="\t")

            with open(pdata_filename, "rb") as file:
                pdict = pickle.load(file)

            num_pts = hf["neutral_ll"].shape[0]
            if cond_only:
                with open(bh_filename, "rb") as file:
                    bf = pickle.load(file)
                idx_list = np.where(bf["bh_classes"] == sel_types.index(sel_type))[0]
                num_pts = idx_list.shape[0]
            elif init_dist == "real_special":
                num_pts = 500
                idx_list = np.arange(num_pts)
            else:
                idx_list = np.arange(num_pts)
            if sel_type == "neutral":
                for run_type in run_types:
                    if run_type != "het":
                        run_type_val = run_type
                    else:
                        run_type_val = "over"

                    if cond_only:
                        if run_type in ["add", "dom", "rec"]:
                            idx_list = np.where(bf["bh_classes"] == sel_types.index(run_type))[0]
                            num_pts = idx_list.shape[0]
                            if run_type == "add" and num_pts == 0:
                                s_vals.extend([10])
                                s_types.extend(["add"])
                                s_strs.extend(["Neutral"])
                                subtract_one_from_first_counts = True
                        else:
                            continue
                    if run_type == "rec" and init_dist in [.005, "recip"]:
                        continue
                    exit_codes = hf[f"{run_type}_run"]["exit_codes"]
                    illegal_s += exit_codes[exit_codes == 2].sum()
                    illegal_s += exit_codes[exit_codes == 4].sum()
                    illegal_s += exit_codes[exit_codes == 12].sum()
                    s_data = get_1d_s_data_from_type(hf[f"{run_type}_run"]["s_final"][:, idx_list[:num_pts]], run_type)
                    if s_data.shape[0] > 0:
                        max_quantile = max(max_quantile, np.quantile(s_data, .99))
                        min_quantile = min(min_quantile, np.quantile(s_data, .01))
                    s_vals.extend(s_data.tolist())
                    s_types.extend([f"{run_type_val}"] * num_pts)
                    s_strs.extend(["Neutral"] * num_pts)
            else:
                if sel_type == "under":
                    s_data = get_1d_s_data_from_type(-hf[f"het_run"]["s_final"][:, idx_list[:num_pts]], sel_type)
                    if s_data.shape[0] > 0:
                        max_quantile = max(max_quantile, np.quantile(s_data, .99))
                        min_quantile = min(min_quantile, np.quantile(s_data, .01))
                    s_vals.extend(s_data.tolist())
                elif sel_type == "over":
                    s_data = get_1d_s_data_from_type(hf[f"het_run"]["s_final"][:, idx_list[:num_pts]], sel_type)
                    if s_data.shape[0] > 0:
                        max_quantile = max(max_quantile, np.quantile(s_data, .99))
                        min_quantile = min(min_quantile, np.quantile(s_data, .01))
                    s_vals.extend(s_data.tolist())
                else:
                    s_data = get_1d_s_data_from_type(hf[f"{sel_type}_run"]["s_final"][:, idx_list[:num_pts]], sel_type)
                    if s_data.shape[0] > 0:
                        max_quantile = max(max_quantile, np.quantile(s_data, .99))
                        min_quantile = min(min_quantile, np.quantile(s_data, .01))
                    s_vals.extend(s_data.tolist())
                s_types.extend([f"{sel_type}"] * num_pts)
                s_strs.extend([sel_str] * num_pts)

        s_types = convert_from_abbrevs(s_types, shortall=True)

        massaged_data = zip(s_vals, s_strs, s_types)

        s_stuff = pd.DataFrame(data=massaged_data, columns=[r"$\hat{s}$", r"True $s$", "Mode of selection"])
        if cond_only:
            box = sns.stripplot(data=s_stuff, x=r"True $s$", y=r"$\hat{s}$", hue="Mode of selection", dodge=True,ax=axs, size=2)
            counts = []
            y_locs = []
            x_locs = []
            for s_i, strength in enumerate(["Neutral", *sel_strs]):
                for m_i, mode in enumerate(["Add.", "Dom.", "Rec."]):
                    bin_mask = ((s_stuff[r"True $s$"]==strength)&(s_stuff["Mode of selection"]==mode)&(s_stuff[r"$\hat{s}$"]<10))
                    counts.append(bin_mask.sum())
                    if counts[-1] > 0:
                        y_locs.append(s_stuff[bin_mask][r"$\hat{s}$"].max())
                    else:
                        y_locs.append(-1)
                    x_locs.append(s_i+(m_i-1)*.28)
                    if counts[-1] > 99:
                        x_locs[-1] += (m_i-1)*.1
            y_locs[0] = y_locs[1]
            for y_i in range(len(y_locs)):
                if y_i%3==2:
                    max_yi = min(max(y_locs[y_i], y_locs[y_i-1], y_locs[y_i-2]), max_quantile)
                    y_locs[y_i] = y_locs[y_i-1] = y_locs[y_i-2] = max_yi*1.05
            for c_i, count in enumerate(counts):
                axs.text(x_locs[c_i], y_locs[c_i], str(count), ha="center", va="bottom", color=colorlist[c_i%3], size=7)
        else:
            box = sns.boxplot(data=s_stuff, x=r"True $s$", y=r"$\hat{s}$", hue="Mode of selection", dodge=True, width=.75,
                  ax=axs,fliersize=1, boxprops={"lw": .5}, medianprops={"lw": .5}, whiskerprops={"lw": .5},
                  capprops={"lw": .5}, flierprops={"alpha": .7}, whis=(2.5, 97.5))
            counts = []
            x_coords = []
            for s_i, strength in enumerate(["Neutral", *sel_strs]):
                for m_i, mode in enumerate(["Add.", "Dom.", "Rec.", "Over.", "Under."]):
                    bin_mask = ((s_stuff[r"True $s$"]==strength)&(s_stuff["Mode of selection"]==mode)&((s_stuff[r"$\hat{s}$"]>max_quantile)|(s_stuff[r"$\hat{s}$"]<min_quantile)))
                    if (bin_counts := bin_mask.sum()) > 0:
                        counts.append(bin_counts)
                    else:
                        counts.append("")
                    x_coords.append(s_i+(m_i-2)*.15)
            for c_i, count in enumerate(counts):
                axs.text(x_coords[c_i], max_quantile*1.01, str(count), ha="center", va="bottom", color=colorlist[c_i%len(colorlist)], size=7)
        axs.axhline(0, color="r", alpha=.65, ls="--", lw=.5)
        for sel_str in sel_strs:
            axs.axhline(sel_str, color="r", alpha=.6, ls="--", lw=.5)
        legend_loc = "lower right" if cond_only else "upper left"
        if not cond_only:
            axs.text(-.24, 1.01, r"$\bf{A}$", fontsize=13, transform=axs.transAxes)
        axs.legend(loc=legend_loc, fontsize=7, labelspacing=.2, handlelength=1.5, handleheight=.5, handletextpad=.4, borderpad=.2, borderaxespad=.2, markerscale=.25 if cond_only else 1)
        axs.set_ylim([min_quantile, max_quantile])
        plt.savefig(f"{output_dir}/g{num_gens}_d{init_dist}_{'strip' if cond_only else 'box'}plots.pdf", format="pdf", bbox_inches="tight")
        plt.close(fig)