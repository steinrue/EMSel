import numpy as np
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
import colorcet as cc

###### MODIFY

data_dir = "data"
output_dir = "output"

###### DO NOT MODIFY

plt.rcParams.update({'font.size': 9,
                     'text.usetex': False,
                     'font.family': 'serif',
                     'font.serif': 'cmr10',
                     'mathtext.fontset': 'cm',
                     'axes.unicode_minus': False,
                     'axes.formatter.use_mathtext': True,})

YEARS_PER_GEN = 8
cmap = cc.m_bkr

binned_path = Path(f"{data_dir}/horse_data_full.csv")
binned_df = np.loadtxt(binned_path, delimiter="\t")
b_st = 22000-binned_df[0, ::3]*YEARS_PER_GEN
b_ns = binned_df[0, 1::3]
b_ff = binned_df[0, 2::3]
fig, axs = plt.subplots(1,1,figsize=(3.1, 3.1), layout="constrained")
axs.invert_xaxis()
zorder=2.1
ls = "-"
lw = 1.7
label = "ASIP"
alpha=1
axs.plot(b_st, b_ff/b_ns, color="k", ls=ls, lw=lw, zorder=zorder, alpha=alpha, label="ASIP")
axs.scatter(b_st, b_ff/b_ns, color="k", marker="o", s=b_ns/1.5, alpha=alpha, zorder=zorder)
for true_size in [10, 25, 40]:
    axs.scatter([], [], s=true_size/1.5, color="k", label=rf"$n_t = {{{true_size}}}$")
axstext = r"$\bf{A}$"
axs.text(-.2, .97, axstext,fontsize=13,transform=axs.transAxes)
axs.legend()
axs.set_ylim([-0.05,1.05])
axs.axhline(0, lw=.5, ls="--", color="k", alpha=.5)
axs.axhline(1, lw=.5, ls="--", color="k", alpha=.5)
axs.set_xlim([1.02*22000, .8*2500])
axs.set_xlabel("Years BP")
axs.set_ylabel("Frequency")
fig.savefig(Path(f"{output_dir}/asip_traj.pdf"), format="pdf", bbox_inches="tight")