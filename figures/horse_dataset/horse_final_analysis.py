import matplotlib.pyplot as plt
from pathlib import Path
import pickle
from cycler import cycler
import pandas as pd

from scipy.stats import chi2
import numpy as np
from emsel.emsel_util import convert_from_abbrevs, bh_correct, get_1d_s_data_from_type, get_llg_array, get_llgka_array, full_bh_procedure, classify_full_run

###### MODIFY

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
colorlist = ["#1D6996", *[coolormap(i) for i in [1,0]], colors[3], colors[4], colors[5]]
plt.rcParams["axes.prop_cycle"] = cycler(color=colorlist)

suffixes = ["full"]
labels = ["B", "C"]
init_ests = ["neut", "uncon"]
k_opt = 1
chi2_dist = chi2(1)

onep_types = ["add", "dom", "rec", "het"]

all_types = ["add", "dom", "rec", "het", "full"]

fancy_col_names = ["Add.\ ($s_2$)", "Dom.\ ($s_2$)", "Rec.\ ($s_2$)", "Het.\ diff. ($s_1$)", "Uncons.\ $(s_1, s_2)$"]

for suffix in suffixes:
    with open(Path(f"{EM_dir}/horse_data_{suffix}_EM.pkl"), "rb") as file:
        hf = pickle.load(file)
    n_ll = hf["neutral_ll"][0]
    for sel_type in onep_types:
        nn_ll = hf[f"{sel_type}_run"]["ll_final"][0]
        llr = 2 * (nn_ll - n_ll)
        p_val = -chi2_dist.logsf(llr) / np.log(10)
        print(f"{suffix} {sel_type} llr: {llr:.3f} logp-value: {p_val:.3f} p-val {np.power(10, -p_val):.8f}")
    full_ll = hf["full_run"]["ll_final"][0]
    d2_stat = 2 * (full_ll - n_ll)
    for i, init_est in enumerate(init_ests):
        with open(Path(f"{EM_dir}/{suffix}_Ne16000_{init_est}_final_sims_EM.pkl"), "rb") as file:
            bootf = pickle.load(file)

        boot_n_lls = bootf["neutral_ll"]
        boot_full_lls = bootf["full_run"]["ll_final"]
        boot_d2s = 2*(boot_full_lls-boot_n_lls)
        print(f"{suffix} empirical percentile: {(boot_d2s<d2_stat).mean():.3f}")
        fig, axs = plt.subplots(1,1,figsize=(3.1, 1.45), layout="constrained")
        axs.hist(boot_d2s, bins=np.arange(np.max(boot_d2s)+1), label="Simulations")
        axs.axvline(d2_stat, color="red", ls="--", label="ASIP $D_2$ value")
        axs.set_xlabel(r"$D_2$ statistic")
        axs.set_ylabel("Counts")
        axs.set_ylim([0, 700])
        axs.legend()
        axs.text(-.22, .95, rf"$\bf{{{labels[i]}}}$", fontsize=13, transform=axs.transAxes)
        fig.savefig(f"{output_dir}/{suffix}_{init_est}_bootstrapped_distribution.pdf", format="pdf", bbox_inches="tight")


with open(Path(f"{EM_dir}/horse_data_full_EM.pkl"), "rb") as file:
    full_hf = pickle.load(file)
full_n_ll = full_hf["neutral_ll"][0]

with open(Path(f"{EM_dir}/horse_data_colon3_EM.pkl"), "rb") as file:
    trunc_hf = pickle.load(file)
trunc_n_ll = trunc_hf["neutral_ll"][0]

full_s_list = []
full_ll_list = []
full_p_list = []
trunc_s_list = []
trunc_ll_list = []
trunc_p_list = []
for final_s_type in all_types:
    full_nn_ll = full_hf[f"{final_s_type}_run"]["ll_final"][0]
    full_llr = (full_nn_ll - full_n_ll)

    trunc_nn_ll = trunc_hf[f"{final_s_type}_run"]["ll_final"][0]
    trunc_llr = (trunc_nn_ll - trunc_n_ll)
    if final_s_type != "full":
        full_p_val = -chi2_dist.logsf(2*full_llr) / np.log(10)
        full_s = get_1d_s_data_from_type(full_hf[f"{final_s_type}_run"]["s_final"], final_s_type)[0]
        trunc_p_val = -chi2_dist.logsf(2*trunc_llr) / np.log(10)
        trunc_s = get_1d_s_data_from_type(trunc_hf[f"{final_s_type}_run"]["s_final"], final_s_type)[0]
    else:
        full_p_val = trunc_p_val = "--"
        full_s_temp = full_hf[f"{final_s_type}_run"]["s_final"][0]
        full_s = f"({full_s_temp[0]:.4f}, {full_s_temp[1]:.4f})"
        trunc_s_temp = trunc_hf[f"{final_s_type}_run"]["s_final"][0]
        trunc_s = f"({trunc_s_temp[0]:.4f}, {trunc_s_temp[1]:.4f})"
        full_llr -= 1
        trunc_llr -= 1
    full_s_list.append(full_s)
    full_ll_list.append(full_llr)
    full_p_list.append(full_p_val)
    trunc_s_list.append(trunc_s)
    trunc_ll_list.append(trunc_llr)
    trunc_p_list.append(trunc_p_val)

full_s_list = [f"{v:.4f}" if isinstance(v, float) else v for v in full_s_list ]
full_ll_list = [f"{v:.2f}" for v in full_ll_list]
full_p_list = [f"{np.power(10, -v):.2e}" if isinstance(v, float) else v for v in full_p_list]
full_p_list = [v.replace("e", " $\cdot\ 10^{").replace("-0", "-") + "}$" if v != "--" else v for v in full_p_list]
trunc_s_list = [f"{v:.4f}" if isinstance(v, float) else v for v in trunc_s_list ]
trunc_ll_list = [f"{v:.2f}" for v in trunc_ll_list]
trunc_p_list = [f"{np.power(10, -v):.2e}" if isinstance(v, float) else v for v in trunc_p_list]
trunc_p_list = [v.replace("e", " $\cdot\ 10^{").replace("-0", "-") + "}$" if v != "--" else v for v in trunc_p_list]
full_ll_list[-1] += "$^\dagger$"
trunc_ll_list[-1] += "$^\dagger$"
asip_table = pd.DataFrame(data=[full_s_list, full_ll_list, full_p_list, trunc_s_list, trunc_ll_list, trunc_p_list], columns = fancy_col_names)

asip_table["  "] = ["\multirow{3}{*}{Full}", "", "", "\multirow{3}{*}{Trunc.}", "", ""]
asip_table[" "] = ["$\hat{s}$", "$\ell\ell - \ell\ell_0$", "p-value"] *2
cols = list(asip_table.columns)
cols = [*cols[-2:]] + cols[:-2]
asip_table = asip_table[cols]
asip_latex = asip_table.to_latex(index=False, column_format="ccccccc")
asip_lines = asip_latex.splitlines()
asip_lines.insert(len(asip_lines)-5, '\midrule')
asip_latex = '\n'.join(asip_lines)
with open(f"{output_dir}/asip_table.tex", "w") as file:
    file.write(asip_latex)




