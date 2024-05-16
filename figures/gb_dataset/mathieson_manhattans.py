import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pickle
from scipy.stats import chi2
from emsel.emsel_util import bh_correct, extendedFisher, windowMatrix

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
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

coolormap = plt.get_cmap("Dark2")
colorlist = ["#1D6996", *[coolormap(i) for i in [1,0]], colors[3], colors[4], colors[5]]


alpha=.05
mathieson_windows = [[157606928, 157803918],[71142350,71183690],[28334927,28526228]]
mathieson_chrs = [7, 11, 15]
axstexts = ["B", "A", "A"]


start_pos = 0
complete_agg_data_path = Path(f"{output_dir}/GB_v54.1_{genodata_type}_agg_data.pkl")
with open(complete_agg_data_path, "rb") as file:
    cdata = pickle.load(file)

all_windows = {}
with open(Path(f"{output_dir}/GB_v54.1_{genodata_type}_add_bh.pkl"), "rb") as file:
    snp_df = pickle.load(file)

windowed_p = windowMatrix(cdata["all_p"][f"add_p"], 25)
brown_p = extendedFisher(windowed_p, standardBrown=True)
brown_p = np.concatenate((np.zeros(26), brown_p, np.zeros(25)))
bp_bh = -np.log10(bh_correct(np.power(10, -brown_p), alpha)[0])


for m_i, mwindow in enumerate(mathieson_windows):
    chrom = mathieson_chrs[m_i]
    l_loc_index, r_loc_index = np.searchsorted(cdata["all_loc_per_chrom"][cdata["all_chrom"]==chrom], [mwindow[0], mwindow[1]])
    l_loc_index += cdata["all_loc"][f"chr_{chrom}_idx_offset"]
    r_loc_index += cdata["all_loc"][f"chr_{chrom}_idx_offset"]
    snp_fig, snp_axs = plt.subplots(1,1,figsize=(3.1,3.1), layout="constrained")
    snp_axs.plot(cdata["all_loc_per_chrom"][l_loc_index-1:r_loc_index+1], cdata["all_p"][f"add_p"][l_loc_index-1:r_loc_index+1], "o", ms=2, c=colorlist[0], label="Raw", zorder=2.1)
    snp_axs.set_ylim([0, 18])
    snp_axs.set_xlabel(f"Position (Mbp) on Chr. {chrom}")
    snp_axs.set_ylabel(r"$-\log_{10}(p)$")
    snp_axs.plot(cdata["all_loc_per_chrom"][l_loc_index-1:r_loc_index+1],
                 brown_p[l_loc_index-1:r_loc_index+1],"^", ms=2, c=colorlist[1], label="Post")
    snp_axs.axhline(-np.log10(snp_df["p_bh"]), color=colorlist[0], ls="--", lw=.75, label=r"Raw BH thresh.")
    snp_axs.axhline(bp_bh, color=colorlist[1], ls="--", lw=.75, label="Post BH thresh.")
    snp_axs.ticklabel_format(axis="x", scilimits=(6, 6))
    plt.setp(snp_axs.get_xaxis().get_offset_text(), visible=False)
    snp_axs.text(-.2, .97, rf"$\bf{{{axstexts[m_i]}}}$", fontsize=13, transform=snp_axs.transAxes)
    snp_axs.legend()
    snp_fig.savefig(f"{output_dir}/{genodata_type}_add_brown_mathieson_c{chrom}_{l_loc_index}.pdf", format="pdf", bbox_inches="tight")
    plt.close(snp_fig)