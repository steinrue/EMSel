import numpy as np
import pickle
from pathlib import Path
from emsel.emsel_util import get_1d_s_data_from_type
import matplotlib.pyplot as plt
###### MODIFY

data_dir = "data/bootstrap"
EM_dir = "EM/bootstrap"
output_dir = "output/bootstrap"
genodata_type = "capture_only"

###### DO NOT MODIFY

c_types = ["add", "dom", "het", "rec"]

for sel_type in c_types:
    s_sets = []
    low_ivl = []
    high_ivl = []
    true_vals = []
    rsids = []
    meds = []
    sds = []
    for fstr in Path(f"{EM_dir}").iterdir():
        if fstr.suffix == ".pkl":
            if sel_type in fstr.stem:
                with open(fstr, "rb") as file:
                    temp_hf = pickle.load(file)
                with open(f"{data_dir}/{fstr.stem[:-3]}_bootstrap_pd.pkl", "rb") as file:
                    temp_pd = pickle.load(file)
                true_vals.append(temp_pd["s2_true"] if sel_type != "het" else temp_pd["s1_true"])
                s_data = get_1d_s_data_from_type(temp_hf[f'{sel_type}_run']["s_final"], sel_type)
                s_sets.append(s_data)
                meds.append(np.median(s_data))
                low_ivl.append(np.quantile(s_data, .025))
                high_ivl.append(np.quantile(s_data, .975))
                rsids.append(fstr.stem.split("_")[0])

    meds = np.array(meds)
    low_ivl = np.array(low_ivl)
    high_ivl = np.array(high_ivl)
    true_vals = np.array(true_vals)
    corr_meds = 2*true_vals-meds
    low_ivl += true_vals-meds
    high_ivl += true_vals-meds

    fig, axs = plt.subplots(1,1,figsize=(2*len(s_sets),5), layout="constrained")
    axs.boxplot(s_sets, whis=(2.5, 97.5))
    axs.plot(np.arange(1,len(s_sets)+1), true_vals, "rs")
    for i in range(len(s_sets)):
        axs.text(i/(len(s_sets)), 1.03, f"{corr_meds[i]:.4f} ({low_ivl[i]:.4f}, {high_ivl[i]:.4f})", transform=axs.transAxes)
        axs.text(i/(len(s_sets)), 1.08, f"{rsids[i]}:", transform=axs.transAxes)

    fig.savefig(Path(f"{output_dir}/{genodata_type}_{sel_type}_bootstrap_boxplots.pdf"), format="pdf", bbox_inches="tight")
    plt.close(fig)
