import numpy as np
import pickle
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

###### MODIFY

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

c_types = ["add", "dom", "het", "rec"]
ticklabels = ["Add.", "Dom.", "Rec.", "Het."]

complete_agg_data_path = Path(f"{output_dir}/GB_v54.1_{genodata_type}_agg_data.pkl")
with open(complete_agg_data_path, "rb") as file:
    cdata = pickle.load(file)

cutoff_val = np.quantile(cdata["all_ll"]["add_ll"], .99)
add_mask = cdata["all_ll"]["add_ll"] > cutoff_val

signed_lls = {}
unsigned_lls = {}

for c_type in c_types:
    signed_lls[f"{c_type}_ll"] = cdata["all_ll"][f"{c_type}_ll"]
    unsigned_lls[f"{c_type}_ll"] = np.copy(signed_lls[f"{c_type}_ll"])
    signed_lls[f"{c_type}_ll"][cdata["all_s"][f"{c_type}_s"] < 0] *= -1
    signed_lls[f"{c_type}_ll"] = signed_lls[f"{c_type}_ll"][add_mask]
sgn_ll_df = pd.DataFrame(signed_lls)
sgn_ll_df = sgn_ll_df[["add_ll", "dom_ll", "rec_ll", "het_ll"]]

sgn_lls_corr = sgn_ll_df.corr(method="spearman")
sgn_lls_corr = sgn_lls_corr[1:-1].drop(columns=["het_ll", "rec_ll"])

sgn_mask = np.tril(np.ones_like(sgn_lls_corr, dtype=bool))
sgn_mask = ~sgn_mask

usgn_ll_df = pd.DataFrame(unsigned_lls)
usgn_ll_df = usgn_ll_df[["add_ll", "dom_ll", "rec_ll", "het_ll"]]

usgn_lls_corr = usgn_ll_df.corr(method="spearman")
usgn_lls_corr = usgn_lls_corr[3:].drop(columns=["het_ll"])

vmin = min(np.amin(usgn_lls_corr.to_numpy()), np.amin(sgn_lls_corr.to_numpy()))
vmax = max(np.amax(usgn_lls_corr.to_numpy()), np.amax(sgn_lls_corr.to_numpy()))

fig, axs = plt.subplots(1,1,figsize=(1.8,1.2), layout="constrained")
sns.heatmap(sgn_lls_corr, cmap="crest_r", linewidth=.8, vmax=vmax, vmin=vmin, cbar=False, fmt=".2f", annot=True, annot_kws={'fontsize':12},mask=sgn_mask,ax=axs, xticklabels=["Add.", "Dom."], yticklabels=["Dom.", "Rec."])
axs.text(-.2, .97, r"$\bf{A}$", fontsize=13, transform=axs.transAxes)
axs.tick_params(axis='both', which='both', length=0, labeltop=False, labelbottom=True)
fig.savefig(f"{output_dir}/{genodata_type}_top_add_signed_ll_correlations.pdf", format="pdf", bbox_inches="tight")

fig2, axs2 = plt.subplots(1,1,figsize=(2.5, .8), layout="constrained")
sns.heatmap(usgn_lls_corr, cmap="crest_r", linewidth=.8,vmax=vmax, vmin=vmin, cbar=False, fmt=".2f", annot=True, annot_kws={'fontsize':12},ax=axs2, xticklabels=["Add.", "Dom.", "Rec."], yticklabels=["Het."])
axs2.text(-.12, .97, r"$\bf{B}$", fontsize=13, transform=axs2.transAxes)
axs2.tick_params(axis='both', which='both', length=0, labeltop=False, labelbottom=True)
fig2.savefig(f"{output_dir}/{genodata_type}_top_add_unsigned_ll_correlations.pdf", format="pdf", bbox_inches="tight")