import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pickle
from emsel.emsel_util import params_dict_to_str, bh_correct, get_roc_and_auc, convert_from_abbrevs, plot_qq
from copy import deepcopy
from scipy.stats import chi2, gengamma
from pandas import DataFrame
import seaborn as sns

###### MODIFY

EM_dir = "EM/pure_sim"
output_dir = "output/pure_sim"

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

subdir_name = "tree/sim"

#chroms = range(1,25)
classification_types = ["neutral", "add", "dom", "rec", "over", "under", "full"]
sel_types = ["add", "dom", "rec", "over", "under"]
sel_types_rows = ["add", "dom", "rec", "over", "under"]
alt_sel_types_rows = ["add", "dom", "rec", "het"]
sel_strs = [.005, .01, .025, .05]

num_gens = 251
init_dist = .25

manhattan_path = Path(f"../runs/{subdir_name}/figures")
auc_df = np.zeros((len(sel_types), len(sel_strs)))

ndict = {}
ndict["sel_type"] = "neutral"
ndict["num_gens"] = num_gens
ndict["init_dist"] = init_dist

neutral_filename = params_dict_to_str(**ndict)
neutral_hmm_path = Path(f"{EM_dir}/{neutral_filename}_EM.pkl")
with open(neutral_hmm_path, "rb") as file:
    nf = pickle.load(file)

neutral_ll = nf["neutral_ll"]
logps = []
labels = []
for type_i, run_type in enumerate(sel_types):
    run_EM_str = run_type if run_type not in ['over', 'under'] else 'het'
    current_bigtable_list = []
    pdict = deepcopy(ndict)
    pdict["sel_type"] = f"{run_type}"
    run_ll = nf[f"add_run"]["ll_final"]
    llr = 2*(run_ll-neutral_ll)
    llr[llr <= 0] = 1e-12
    llr_med = np.median(llr)

    gengamma_sl_fit = gengamma(*gengamma.fit(llr[llr > llr_med] - llr_med, floc=0, fscale=1))
    chisq_sl_fit = chi2(1)
    full_p_vals = -chisq_sl_fit.logsf(llr)/np.log(10)

    for str_i, sel_str in enumerate(sel_strs):
        pdict["sel_str"] = sel_str
        nn_fname = params_dict_to_str(**pdict)
        nn_path = Path(f"{EM_dir}/{nn_fname}_EM.pkl")
        if not nn_path.is_file():
            print(f"{nn_path} not found")
            sf_llr = np.zeros(nf["neutral_ll"].shape[0])
        else:
            with open(nn_path, "rb") as file:
                sf = pickle.load(file)
            sf_run_ll = sf[f"add_run"]["ll_final"]
            sf_llr = 2*(sf_run_ll-sf["neutral_ll"])
            print(f"{nn_fname} {run_type} stats: {np.min(sf_llr):.4f} {(sf_llr < 0).sum()}")
        nn_llr = sf_llr
        nn_p_vals = -chisq_sl_fit.logsf(nn_llr)/np.log(10)
        roc_FPR, roc_TPR, auc = get_roc_and_auc(np.power(10, -full_p_vals), np.power(10, -nn_p_vals))
        auc_df[type_i, str_i] = f"{auc:.2f}"
    if run_type == "under":
        continue
    logps.append(full_p_vals)
    labels.append(convert_from_abbrevs(run_EM_str))

big_ts_dict = {}
big_ts_dict["num_gens"] = f"{num_gens}"
big_ts_dict["init_dist"] = f"{init_dist}"
big_ts_str = params_dict_to_str(**big_ts_dict)
big_df = DataFrame(auc_df, columns=sel_strs, index=convert_from_abbrevs(sel_types_rows))
big_df = big_df.loc[~(big_df == 0).all(axis=1)]
fig_height = 4.5 if big_df.shape[0] == len(sel_types_rows) else 3.75
fig_width = 3.1
add_em_auc_fig, add_em_auc_axs = plt.subplots(1, 1, figsize=(fig_width, fig_height/(8/fig_width)), layout="constrained")
add_em_auc_axs.text(-.55, 1.05, r"$\bf{C}$", fontsize=13, transform=add_em_auc_axs.transAxes)
sns.heatmap(big_df, cmap="crest_r", linewidth=.8, cbar=False, fmt=".2f", annot=True, annot_kws={'fontsize':10}, ax=add_em_auc_axs)
add_em_auc_axs.tick_params(axis='both', which='both', length=0, labeltop=True, labelbottom=False)
add_em_auc_axs.set_xlabel(r"$\bf{s}$", fontsize=10)
add_em_auc_axs.xaxis.set_label_position('top')
add_em_auc_axs.set_ylabel(r"$\bf{Mode}$", fontsize=10)
add_em_auc_fig.savefig(Path(f"{output_dir}/{big_ts_str}_add_em_auc_plot.pdf"), format="pdf", bbox_inches="tight")
plt.close(add_em_auc_fig)

auc_df = np.zeros((len(sel_types), len(sel_strs)))
logps = []
labels = []
for type_i, run_type in enumerate(sel_types):
    run_EM_str = run_type if run_type not in ['over', 'under'] else 'het'
    if run_type == "under":
        auc_df = auc_df[:-1, :]
        continue
    current_bigtable_list = []
    pdict = deepcopy(ndict)
    pdict["sel_type"] = "add"
    run_ll = nf[f"{run_EM_str}_run"]["ll_final"]
    llr = 2*(run_ll-neutral_ll)
    llr[llr <= 0] = 1e-12
    llr_med = np.median(llr)

    gengamma_sl_fit = gengamma(*gengamma.fit(llr[llr > llr_med] - llr_med, floc=0, fscale=1))
    chisq_sl_fit = chi2(1)

    med_p_vals = np.zeros_like(llr)
    med_p_vals[llr > llr_med] = (1 - gengamma_sl_fit.cdf(llr[llr > llr_med] - llr_med)) / 2
    med_p_vals[llr <= llr_med] = np.clip(1 - llr[llr <= llr_med] / (2 * llr_med), .5, 1)

    full_p_vals = -chisq_sl_fit.logsf(llr)/np.log(10)

    for str_i, sel_str in enumerate(sel_strs):
        pdict["sel_str"] = sel_str
        nn_fname = params_dict_to_str(**pdict)
        nn_path = Path(f"{EM_dir}/{nn_fname}_EM.pkl")
        if not nn_path.is_file():
            print(f"{nn_path} not found")
            sf_llr = np.zeros(nf["neutral_ll"].shape[0])
        else:
            with open(nn_path, "rb") as file:
                sf = pickle.load(file)
            sf_run_ll = sf[f"{run_EM_str}_run"]["ll_final"]
            sf_llr = 2*(sf_run_ll-sf["neutral_ll"])
            print(f"{nn_fname} {run_type} stats: {np.min(sf_llr):.4f} {(sf_llr < 0).sum()}")
        nn_llr = sf_llr
        nn_p_vals = -chisq_sl_fit.logsf(nn_llr)/np.log(10)
        roc_FPR, roc_TPR, auc = get_roc_and_auc(np.power(10, -full_p_vals), np.power(10, -nn_p_vals))
        auc_df[type_i, str_i] = f"{auc:.2f}"
    logps.append(full_p_vals)
    labels.append(convert_from_abbrevs(run_EM_str))

big_ts_dict = {}
big_ts_dict["num_gens"] = f"{num_gens}"
big_ts_dict["init_dist"] = f"{init_dist}"
big_ts_str = params_dict_to_str(**big_ts_dict)
big_df = DataFrame(auc_df, columns=sel_strs, index=convert_from_abbrevs(alt_sel_types_rows, shorthet=True))
big_df = big_df.loc[~(big_df == 0).all(axis=1)]
fig_height = 4.5 if big_df.shape[0] == len(sel_types_rows) else 3.75
fig_width = 3.1
add_sim_auc_fig, add_sim_auc_axs = plt.subplots(1, 1, figsize=(fig_width, fig_height/(8/fig_width)), layout="constrained")
sns.heatmap(big_df, cmap="crest_r", linewidth=.8, cbar=False, fmt=".2f", annot=True, annot_kws={'fontsize':10}, ax=add_sim_auc_axs)
add_sim_auc_axs.text(-.35, 1.2, r"$\bf{D}$", fontsize=13, transform=add_sim_auc_axs.transAxes)
add_sim_auc_axs.tick_params(axis='both', which='both', length=0, labeltop=True, labelbottom=False)
add_sim_auc_axs.set_xlabel(r"$\bf{s}$", fontsize=10)
add_sim_auc_axs.xaxis.set_label_position('top')
add_sim_auc_axs.set_ylabel(r"$\bf{Mode}$", fontsize=10)
add_sim_auc_fig.savefig(Path(f"{output_dir}/{big_ts_str}_add_sim_auc_plot.pdf"), format="pdf", bbox_inches="tight")
plt.close(add_sim_auc_fig)