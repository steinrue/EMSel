import numpy as np
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
import colorcet as cc

###### MODIFY

data_dir = "data"
EM_dir = "EM"
output_dir = "output"
genodata_type = "capture_only"
classification_types = ["add", "dom", "rec", "het"]

###### DO NOT MODIFY

plt.rcParams.update({'font.size': 9,
                     'text.usetex': False,
                     'font.family': 'serif',
                     'font.serif': 'cmr10',
                     'mathtext.fontset': 'cm',
                     'axes.unicode_minus': False,
                     'axes.formatter.use_mathtext': True,})

#these are the lead SNPs for the non-significant example on chromosome 7 and the comparison to Mathieson and Terhorst
lead_snp_idxs = np.array([355758, 490617, 597996])
region_size = 20
YEARS_PER_GEN = 28.1

cmap = cc.m_bkr

complete_agg_data_path = Path(f"{output_dir}/GB_v54.1_{genodata_type}_agg_data.pkl")
with open(complete_agg_data_path, "rb") as file:
    cdata = pickle.load(file)

for c_type in classification_types:
    with open(f"{output_dir}/{genodata_type}_{c_type}_sig_windows.pkl",
              "rb") as file:
        temp_windows = pickle.load(file)
        for _, tr in temp_windows.iterrows():
            if min(np.abs(int(tr["SNP index of max."])-lead_snp_idxs)) > region_size:
                lead_snp_idxs = np.concatenate((lead_snp_idxs, np.array([int(tr["SNP index of max."])])))
lead_snp_idxs = [[lead_idx] for lead_idx in lead_snp_idxs]

#additional manual additions to HYDIN
lead_snp_idxs[-1] = [*lead_snp_idxs[-1], 635227, 635233]

binned_path = Path(f"{data_dir}/GB_v54.1_{genodata_type}_complete_data_binned.csv")
binned_df = np.loadtxt(binned_path, delimiter="\t")
b_st = 4450-binned_df[:, ::3]*YEARS_PER_GEN
b_ns = binned_df[:, 1::3]
b_ff = binned_df[:, 2::3]

for l_i, lead_snp in enumerate(lead_snp_idxs):
    fig, axs = plt.subplots(1,1,figsize=(3.1, 3.1), layout="constrained")
    axs.invert_xaxis()

    for c_i, snp_idx in enumerate(np.arange(lead_snp[0]-region_size, lead_snp[0]+region_size+1)):
        if snp_idx in lead_snp:
            zorder=2.1
            ls = "-"
            lw = 2
            label = cdata['all_rsid'][snp_idx]
            alpha=1
        else:
            zorder = 2
            ls = "--"
            lw = 1
            label = ""
            alpha=.5
        axs.plot(b_st[snp_idx], b_ff[snp_idx]/b_ns[snp_idx], color=cmap(c_i/(region_size*2+1)), ls=ls, lw=lw, zorder=zorder, alpha=alpha, label=label)
        axs.scatter(b_st[snp_idx], b_ff[snp_idx]/b_ns[snp_idx], color=cmap(c_i/(region_size*2+1)), marker="o", s=b_ns[snp_idx]*.5, alpha=alpha, zorder=zorder)
    if 355758 in lead_snp:
        axstext = r"$\bf{D}$"
    elif 250327 in lead_snp:
        axstext = r"$\bf{C}$"
    else:
        axstext = r"$\bf{B}$"
    axs.text(-.2, .97, axstext,fontsize=13,transform=axs.transAxes)
    axs.legend()
    axs.set_ylim([-0.05,1.05])
    axs.axhline(0, lw=.5, ls="--", color="k", alpha=.5)
    axs.axhline(1, lw=.5, ls="--", color="k", alpha=.5)
    axs.set_xlim(right=999.99)
    axs.set_xlabel("Years BP")
    axs.set_ylabel("Frequency")
    fig.savefig(Path(f"{output_dir}/{genodata_type}_chr{int(cdata['all_chrom'][lead_snp[0]])}_{lead_snp[0]}_{cdata['all_rsid'][lead_snp[0]]}.pdf"), format="pdf", bbox_inches="tight")