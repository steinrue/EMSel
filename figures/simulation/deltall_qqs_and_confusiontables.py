import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from pathlib import Path
import pickle
from emsel_util import classify_full_run, get_one_count_matrix, params_dict_to_str, get_llg_array, get_lr_statistic, get_llgka_array, full_bh_procedure, convert_from_abbrevs, convert_to_abbrevs, plot_qq
from copy import deepcopy
from scipy.stats import chi2, gengamma
from pandas import DataFrame
from scipy.signal import savgol_filter
import seaborn as sns
from cycler import cycler


###### MODIFY
sel_strs = [.005, .01, .025, .05]
num_gens_list = [101, 251, 1001]
init_dists = [.005, .25, "recip"]

save_gengamma = False

EM_dir = "EM"
classified_dir = "classified"
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
colorlist = ["#1D6996", *[coolormap(i) for i in [1,0]], colors[3], colors[4], colors[5]]
init_colorlist = colorlist
plt.rcParams["axes.prop_cycle"] = cycler(color=colorlist)

classification_types = ["neutral", "add", "dom", "rec", "over", "under", "full"]
classification_types_rows = ["Neut.", "Add.", "Dom.", "Rec.", "Over.", "Under.", "Indet."]
sel_types = ["add", "dom", "rec", "over", "under"]
sel_types_rows = ["add", "dom", "rec", "over", "under"]



onep_types = ["add", "dom", "rec"]
update_types = [*onep_types, "full"]
final_sel_types = ["add", "dom", "rec", "over", "under"]

gengamma_save_path = Path(f"{output_dir}/gengamma_params.pkl")

p_val = .05
save_bh = True

bigtable_idxs = []
full_bigtable_list = []
for num_gens in num_gens_list:
    for init_dist in init_dists:
        row_list = []
        g_list = []
        ndict = {}
        ndict["sel_type"] = "neutral"
        ndict["num_gens"] = num_gens
        ndict["init_dist"] = init_dist

        neutral_filename = params_dict_to_str(**ndict)
        row_list.append([neutral_filename])
        neutral_hmm_path = Path(f"{EM_dir}/{neutral_filename}_EM.pkl")
        with open(neutral_hmm_path, "rb") as file:
            nf = pickle.load(file)

        n_ll_array = get_llg_array(nf, onep_types, classify_full_run(nf["full_run"]["s_final"])[0])[:, :, np.newaxis]
        final_sel_types = ["add", "dom", "rec", "over", "under"]
        for sel_str in sel_strs:
            pdict = deepcopy(ndict)
            llg_array_list = []
            row_strs = []
            for sel_i, sel_type in enumerate(sel_types):
                pdict["sel_type"] = sel_type
                pdict["sel_str"] = sel_str
                exp_name = params_dict_to_str(**pdict)
                onep_path = Path(f"{exp_name}_EM.pkl")
                if sel_type == "under" and not onep_path.is_file():
                    final_sel_types = ["add", "dom", "rec", "over"]
                    continue

                with open(onep_path, "rb") as file:
                    ohf = pickle.load(file)

                classified_array, _ = classify_full_run(ohf["full_run"]["s_final"])
                temp_llg_array = get_llg_array(ohf, onep_types, classified_array)
                llg_array_list.append(temp_llg_array)
                row_strs.append(convert_from_abbrevs(sel_type))
            row_list.append(row_strs)
            g_ll_array = np.dstack(llg_array_list)
            g_list.append(g_ll_array)

        usable_length = len(final_sel_types)
        k_opt_array = np.zeros(5000)
        num_opt = 0
        k_opt = 0
        k_range = np.linspace(-5, 5, 5000)
        for i, k in enumerate(k_range):
            final_count_array = np.zeros((len(final_sel_types), len(onep_types) + 3))
            for g_i, g_array in enumerate(g_list):
                temp_count_matrix = get_one_count_matrix(k, g_array)
                final_count_array[:usable_length, :] += temp_count_matrix[:usable_length, :]
            num_temp = np.trace(final_count_array)
            k_opt_array[i] = num_temp
        k_opt_sg = savgol_filter(k_opt_array, 101, 3)
        k_opt = k_range[np.argmax(k_opt_sg)]
        k_opt_noisy = k_range[np.argmax(k_opt_array)]


        #k_opt = 1.5
        nk_array = get_llgka_array(n_ll_array, k_opt, 0)
        n_lrs = get_lr_statistic(nk_array).flatten()
        n_lrs[n_lrs < 0] = 1e-14
        lr_shift = np.median(n_lrs)
        gengamma_sl_fit = gengamma(*gengamma.fit(n_lrs[n_lrs > lr_shift] - lr_shift, floc=0, fscale=1))
        chisq1_sl_fit = chi2(1)
        chisq2_sl_fit = chi2(2)

        if save_gengamma:
            gengamma_dict = {"k_opt": k_opt, "lr_shift": lr_shift, "gengamma_fit": gengamma_sl_fit}
            with open(gengamma_save_path, "wb") as file:
                pickle.dump(gengamma_dict, file)

        med_p_vals = np.zeros_like(n_lrs)
        med_p_vals[n_lrs > lr_shift] = (1 - gengamma_sl_fit.cdf(n_lrs[n_lrs > lr_shift] - lr_shift)) / 2
        med_p_vals[n_lrs <= lr_shift] = np.clip(1 - n_lrs[n_lrs <= lr_shift] / (2 * lr_shift), .5, 1)

        full_p_vals_1 = 1 - chisq1_sl_fit.cdf(n_lrs)
        full_p_vals_2 = 1 - chisq2_sl_fit.cdf(n_lrs)

        blank_array = np.empty(len(classification_types), dtype=str)
        blank_array.fill("")

        gllgka_list = []
        for gla in g_list:
            gllgka_list.append(get_llgka_array(gla, k_opt, 0))
        p_cutoff, g_p_val_matrix, test_vals = full_bh_procedure([nk_array, *gllgka_list], gengamma_sl_fit, lr_shift,
                                                                p_val, bh=False)

        n_full = np.array([np.sum(test_vals[0] == i)/test_vals[0].shape[0] for i in np.arange(len(classification_types))])
        n_classified_vals = [np.sum(test_vals[0] == i)/test_vals[0].shape[0] for i in np.arange(len(classification_types))]
        n_max_idx = np.argmax(n_classified_vals)
        n_full[n_max_idx] = n_full[n_max_idx]
        n_full[0] = n_full[0]

        indices = [f"Neutral"]
        if save_bh:
            tpath = Path(f"{row_list[0][0]}_classified.pkl")
            bh_dict = {
                "bh_classes": test_vals[0],
                "p_vals": g_p_val_matrix[0],
                "p_cutoff": p_cutoff
            }
            with open(tpath, "wb") as file:
                pickle.dump(bh_dict, file)
        test_vals = test_vals[1:]
        g_p_val_matrix = g_p_val_matrix[1:]
        row_list = row_list[1:]
        full_excel_array = np.copy(n_full)
        cut_idx = 0
        params_dict_to_str(**{"init_dist": init_dist, "num_gens": num_gens})
        for s_i, s_array in enumerate(test_vals):
            for col in np.arange(s_array.shape[1]):
                tpath = Path(f"{classified_dir}/{convert_to_abbrevs(row_list[s_i][col])}_s{str(sel_strs[s_i])[2:]}_g{num_gens}_d{str(init_dist)[2:]}_classified.pkl")
                bh_dict = {
                    "bh_classes": s_array[:, col],
                    "p_vals": g_p_val_matrix[s_i][:, col],
                    "p_cutoff": p_cutoff,
                }
                with open(tpath, "wb") as file:
                    pickle.dump(bh_dict, file)
        for s_i, s_array in enumerate(test_vals):
            if sel_strs[s_i] in [0.005, 0.025]:
                continue
            for col in np.arange(s_array.shape[1]):
                t_full_array = np.array(
                    [np.sum(s_array[:, col] == i)/s_array.shape[0] for i in np.arange(len(classification_types))])
                classified_vals = np.array([np.sum(s_array[:, col] == i)/s_array.shape[0] for i in np.arange(len(classification_types))])
                max_idx = np.argmax(classified_vals)
                t_full_array[max_idx] = t_full_array[max_idx]
                t_full_array[col+1] = t_full_array[col+1]
                indices.append(row_list[s_i][col])
                full_excel_array = np.vstack((full_excel_array, t_full_array))

            for row_idx in indices:
                if row_idx == "":
                    continue

        n_full_df = DataFrame(full_excel_array[cut_idx:, :], index=indices[cut_idx:],
                              columns=classification_types_rows)
        n_full_df = n_full_df.loc[~(n_full_df == 0).all(axis=1)]

        vmax = n_full_df.to_numpy().max()
        vmin = n_full_df.to_numpy().min()
        fig2, axs = plt.subplots(3, 1, figsize=(4, 4.25), layout="constrained",
                                 gridspec_kw={'height_ratios': [1, len(final_sel_types), len(final_sel_types)],
                                              'hspace': .1})

        sns.heatmap(n_full_df.iloc[:1, :], vmin=vmin, vmax=vmax, cmap="crest_r", linewidth=.8, cbar=False, fmt=".2f",
                    annot=True, annot_kws={'fontsize': 10}, ax=axs[0])
        sns.heatmap(n_full_df.iloc[1:len(final_sel_types) + 1, :], vmin=vmin, vmax=vmax, cmap="crest_r", linewidth=.8,
                    cbar=False, fmt=".2f", annot=True, annot_kws={'fontsize': 10}, ax=axs[1])
        sns.heatmap(n_full_df.iloc[len(final_sel_types) + 1:, :], vmin=vmin, vmax=vmax, cmap="crest_r", linewidth=.8,
                    cbar=False, fmt=".2f", annot=True, annot_kws={'fontsize': 10}, ax=axs[2])
        axs[0].tick_params(axis='both', which='both', length=0, labeltop=True, labelbottom=False)
        axs[1].tick_params(axis='both', which='both', length=0, labeltop=False, labelbottom=False)
        axs[2].tick_params(axis='both', which='both', length=0, labeltop=False, labelbottom=False)
        axs[0].set_yticklabels(axs[0].get_yticklabels(), rotation=0)
        axs[0].set_xlabel(r"$\bf{Classified}$ $\bf{mode}$", fontsize=10)
        axs[0].xaxis.set_label_position('top')
        axs[1].set_ylabel("s = 0.01")
        axs[2].set_ylabel("s = 0.05")
        fig2.supylabel(r"$\bf{Simulated}$ $\bf{mode}$", fontsize=10)

        offset = 0.05
        edge_col = "orange"
        alpha = .8
        lw = 2
        axs[0].add_patch(
            Rectangle((0 + offset, 0 + offset), 1 - 2 * offset, 1 - 2 * offset, ec=edge_col, fc='none', lw=lw,
                      alpha=alpha))
        for diag_pos in range(len(final_sel_types)):
            axs[1].add_patch(
                Rectangle((diag_pos + 1 + offset, diag_pos + offset), 1 - 2 * offset, 1 - 2 * offset, ec=edge_col,
                          fc='none', lw=lw, alpha=alpha))
            axs[2].add_patch(
                Rectangle((diag_pos + 1 + offset, diag_pos + offset), 1 - 2 * offset, 1 - 2 * offset, ec=edge_col,
                          fc='none', lw=lw, alpha=alpha))
        tablename = params_dict_to_str(**{"init_dist": init_dist, "num_gens": num_gens})
        fig2.savefig(Path(f"{output_dir}/{tablename}_confusion_plot.pdf"),format="pdf", bbox_inches="tight")

        fig, axs = plt.subplots(1, 1, figsize=(3.1, 3.1), layout="constrained")
        logps = [-np.log10(med_p_vals), -np.log10(full_p_vals_1), -np.log10(full_p_vals_2)]
        labels = ["Gen. gamma", r"$\chi^2(1)$", r"$\chi^2(2)$"]
        colors = [colorlist[2], "#861A5E", colorlist[0]]
        axins = axs.inset_axes([.65, .11, .3, .3])
        axs.text(-.2, .97, r"$\bf{D}$", fontsize=13, transform=axs.transAxes)
        plot_qq(axs, axins, logps, labels, colors=colors, legend_loc = "upper left", thin=True)
        fig.savefig(Path(f"{output_dir}/neutral_g{num_gens}_d{init_dist}_dll_all.pdf"), format="pdf", bbox_inches="tight")
        plt.close(fig)