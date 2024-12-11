import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pickle
from emsel.emsel_util import average_p_vals, params_dict_to_str, plot_qq
from scipy.stats import chi2

###### MODIFY

true_EM_dir = "EM/real_matched"
perm_EM_dir = "EM/real_matched/permutations"
output_dir = "output/real_matched"

###### DO NOT MODIFY

Ne = 9987

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

num_gens = 125
init_dist = "real_special"

ndict = {}
ndict["sel_type"] = "neutral"
ndict["num_gens"] = num_gens
ndict["init_dist"] = init_dist
neutral_filename = params_dict_to_str(**ndict)


fig, axs = plt.subplots(1,1,figsize=(3.1,3.1), layout="constrained")
axs.text(-.2, .97, r"$\bf{A}$", fontsize=13, transform=axs.transAxes)
axins = axs.inset_axes([.67, .11, .28, .28])
logps = []
labels = []
for perm_i in range(100):
    neutral_hmm_path = Path(f"{perm_EM_dir}/{neutral_filename}_perm{perm_i}_Ne{Ne}_EM.pkl")
    with open(neutral_hmm_path, "rb") as file:
        nf = pickle.load(file)
    neutral_ll = nf["neutral_ll"]
    run_ll = nf[f"add_run"]["ll_final"]
    llr = 2 * (run_ll - neutral_ll)
    llr[llr <= 0] = 1e-12
    full_p_vals = -chi2(1).logsf(llr) / np.log(10)
    logps.append(full_p_vals)
    labels.append("")

true_neutral_path = Path(f"{true_EM_dir}/{neutral_filename}_Ne{Ne}_EM.pkl")
with open(true_neutral_path, "rb") as file:
    nf = pickle.load(file)
neutral_ll = nf["neutral_ll"]
run_ll = nf[f"add_run"]["ll_final"]
llr = 2 * (run_ll - neutral_ll)
llr[llr <= 0] = 1e-12
full_p_vals = -chi2(1).logsf(llr) / np.log(10)
len_ps = full_p_vals.shape[0]

plot_qq(axs, axins, logps, labels, thin=True, rasterized=True)

axins.plot(np.arange(1, len_ps + 1) / len_ps,
                       np.power(10, -np.sort(full_p_vals)[::-1]), lw=2.25, label="Unpermuted", color="k")
axs.plot(-np.log10(np.arange(1, len_ps + 1) / len_ps),
         np.sort(full_p_vals)[::-1], ".", color="k", markersize=10)
handles, labels = axins.get_legend_handles_labels()
axs.legend(handles, labels, loc="upper right")

fig.savefig(Path(f"{output_dir}/rf_permutations_rasterized.pdf"), format="pdf", bbox_inches="tight", dpi=1000)
plt.close(fig)