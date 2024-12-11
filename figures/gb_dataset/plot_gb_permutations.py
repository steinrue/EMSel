import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pickle
from scipy.stats import chi2
from tqdm import tqdm
from emsel.emsel_util import plot_qq

###### MODIFY

EM_dir = "EM"
output_dir = "output"
genodata_type = "capture_only"

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

chroms = range(1,23)

total_snps = 0
#chrom = 1

all_p_vals = np.array([0])
for chrom in tqdm(chroms):

    em_path = Path(f"{EM_dir}/GB_v54.1_{genodata_type}_c{chrom}_permuted_EM.pkl")
    with open(em_path, "rb") as file:
        hf = pickle.load(file)
    neutral_ll = hf["neutral_ll"]
    run_ll = hf[f"add_run"]["ll_final"]
    llr = 2 * (run_ll - neutral_ll)
    llr[llr <= 0] = 1e-12
    p_vals = -chi2(1).logsf(llr)/np.log(10)
    all_p_vals = np.concatenate((all_p_vals, p_vals))
all_p_vals = all_p_vals[1:]

complete_agg_data_path = Path(f"{output_dir}/GB_v54.1_{genodata_type}_agg_data.pkl")
with open(complete_agg_data_path, "rb") as file:
    agg_data = pickle.load(file)

assert agg_data["all_p"]["add_p"].shape[0] == all_p_vals.shape[0]

fig, axs = plt.subplots(1, 1, figsize=(3.1, 3.1), layout="constrained")
axs.text(-.2, .97, r"$\bf{B}$", fontsize=13, transform=axs.transAxes)
axins = axs.inset_axes([.67, .11, .28, .28])
plot_qq(axs, axins, [agg_data["all_p"][f"add_p"], all_p_vals], ["Unpermuted", "Permuted"], thin=True, rasterized=True)
fig.savefig(f"{output_dir}/{genodata_type}_permuted_qqs_rasterized.pdf", format="pdf", bbox_inches="tight", dpi=600)
plt.close(fig)
