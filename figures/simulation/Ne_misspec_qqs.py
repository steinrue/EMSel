import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pickle
from emsel_util import plot_qq
from scipy.stats import chi2

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

def cmap(val):
    return (1-val, 0, val, 1)
coolormap = plt.get_cmap("Dark2")
colorlist = ["#1D6996", *[coolormap(i) for i in [1,0]], colors[3], colors[4], colors[5]]

Ne_labels = ["5", "10", "20"]
Ne_names = ["Ne5000_", "", "Ne20000_"]
logps = []
labels = []

for Ne_i, Ne_name in enumerate(Ne_names):
    neutral_path = Path(f"{EM_dir}/neutral_g251_d25_{Ne_name}EM.pkl")
    with open(neutral_path, "rb") as file:
        nf = pickle.load(file)

    llr = 2*(nf["add_run"]["ll_final"]-nf["neutral_ll"])
    ps = -chi2(1).logsf(llr)/np.log(10)
    logps.append(ps)
    labels.append(f"$N_e = {Ne_labels[Ne_i]},000$")

fig, axs = plt.subplots(1,1,figsize=(3.1,3.1), layout="constrained")
axins = axs.inset_axes([.67, .11, .28, .28])
plot_qq(axs, axins, logps, labels, legend_loc="upper left")
fig.savefig(Path(f"{output_dir}/Ne_qqs.pdf"), format="pdf", bbox_inches="tight")