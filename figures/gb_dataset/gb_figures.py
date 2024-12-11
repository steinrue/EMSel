import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pickle
from emsel.emsel_util import bh_correct, extendedFisher, windowMatrix, convert_from_abbrevs, plot_qq
import pandas as pd

###### MODIFY

data_dir = "data"
EM_dir = "EM"
output_dir = "output"
genodata_type = "capture_only"
classification_types = ["add", "dom", "rec", "het"]
suffix = ""

###### DO NOT MODIFY

def combine_pos(row):
    return f"{int(row['Left pos.']):,d}-{int(row['Right pos.']):,d}"

def special_format(flt):
    if flt > 1:
        return f"{flt:.2f}"
    return f"{flt:.3f}"

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
colorlist = ["#1d6996", *[coolormap(i) for i in [1,0]], colors[3], colors[4], colors[5]]
full_labels = ["Additive", "Dominant", "Recessive", "Overdom.", "Underdom.", "Indet."]
chroms = range(1,23)
alpha = .05
scatter_markersize = 2
scatter_bh_width = .75
complete_agg_data_path = Path(f"{output_dir}/GB_v54.1_{genodata_type}_{suffix}agg_data.pkl")
onep_classification_types = ["add", "dom", "rec", "het"]
all_classification_types = ["add", "dom", "rec", "het", "full"]
with open(complete_agg_data_path, "rb") as file:
    cdata = pickle.load(file)
all_windows = {}
for c_type in classification_types:
    with open(Path(f"{output_dir}/GB_v54.1_{genodata_type}_{c_type}_bh.pkl"), "rb") as file:
        snp_df = pickle.load(file)
    windowed_p = windowMatrix(cdata["all_p"][f"{c_type}_p"], 25)
    brown_p = extendedFisher(windowed_p, standardBrown=True)
    brown_p = np.concatenate((np.zeros(26), brown_p, np.zeros(25)))
    bp_bh = -np.log10(bh_correct(np.power(10, -brown_p), alpha)[0])
    brown_p_sig_idx = np.where(brown_p > bp_bh)[0]
    brown_diffs = np.diff(brown_p_sig_idx)
    brown_window_boundaries = np.concatenate(([0], np.where(brown_diffs>1)[0]))
    brown_p_windows = []
    for i in range(brown_window_boundaries.shape[0] - 1):
        brown_p_windows.append(
            brown_p_sig_idx[brown_window_boundaries[i] + 1:brown_window_boundaries[i + 1] + 1].tolist())
    brown_p_windows.append(brown_p_sig_idx[brown_window_boundaries[-1] + 1:].tolist())
    valid_windows = []
    for brown_window in brown_p_windows:
        if np.intersect1d(brown_window, snp_df["snp_idx"]).shape[0] > 0:
            if cdata["all_chrom"][brown_window[0]] == cdata["all_chrom"][brown_window[-1]]:
                valid_windows.append(brown_window)
    for v_window in valid_windows:
        chrom = int(cdata["all_chrom"][v_window[0]])
        for zoomscatter_type in ["brown", "plain"]:

            #remove this to generate scatterplots for additional modes of selection - they do not appear in the manuscript or supplement however
            if not ((c_type == "add" and zoomscatter_type == "brown") or (c_type == "full" and zoomscatter_type =="plain")):
                continue
            snp_fig, snp_axs = plt.subplots(1,1,figsize=(3.1, 3.1),layout="constrained")
            if "full" in c_type and "plain" in zoomscatter_type:
                snp_axs.plot(cdata["all_loc_per_chrom"][v_window[0]-2:v_window[-1]+2], cdata["all_p"][f"{c_type}_p"][v_window[0]-2:v_window[-1]+2], "o", markersize=1.5, color=colorlist[0], alpha=.5)
            else:
                snp_axs.plot(cdata["all_loc_per_chrom"][v_window[0] - 2:v_window[-1] + 2],
                             cdata["all_p"][f"{c_type}_p"][v_window[0] - 2:v_window[-1] + 2], "o", ms=scatter_markersize, color=colorlist[0], label="Raw", zorder=2.1)
            snp_axs.set_xlabel(f"Position (Mbp) on Chr. {chrom}")
            snp_axs.set_ylabel(r"$-\log_{10}(p)$")
            if "brown" in zoomscatter_type:
                snp_axs.plot(cdata["all_loc_per_chrom"][v_window[0]-2:v_window[-1]+2],
                             brown_p[v_window[0]-2:v_window[-1]+2],"^", ms=scatter_markersize, c=colorlist[1], label="Post")
            elif "full" in c_type:
                window_mask = (snp_df["snp_idx"]>=v_window[0])&(snp_df["snp_idx"]<=v_window[-1])
                max_loc = np.argmax(snp_df["p"][window_mask])
                max_class = snp_df["classes"][window_mask][max_loc]
                for c_i in range(1,len(full_labels)+1):
                    color_mask = window_mask&(snp_df["classes"]==c_i)
                    if (type_sum := color_mask.sum())>0:
                        snp_axs.plot(snp_df["snp_pos"][color_mask], snp_df["p"][color_mask], "o", ms=scatter_markersize, color=colorlist[c_i-1], label=f"{full_labels[c_i-1]} ({type_sum})")
                    if c_i == max_class:
                        snp_axs.plot(snp_df["snp_pos"][window_mask][max_loc], snp_df["p"][window_mask][max_loc], "*", ms=scatter_markersize*2.5, color=colorlist[c_i-1])
            snp_axs.axhline(-np.log10(snp_df["p_bh"]), color=colorlist[0], ls="--", lw=scatter_bh_width, label="Raw BH thresh.")
            if "brown" in zoomscatter_type:
                snp_axs.axhline(bp_bh, color=colorlist[1], ls="--", lw=scatter_bh_width, label="Post BH thresh.")
            if zoomscatter_type=="plain" and c_type=="full":
                snp_axs.legend(fontsize=7, labelspacing=.2, handlelength=1.5, handleheight=.5,
                           handletextpad=.4, borderpad=.2, borderaxespad=.2)
            else:
                snp_axs.legend()
            snp_axs.ticklabel_format(axis="x",scilimits=(6,6))
            plt.setp(snp_axs.get_xaxis().get_offset_text(), visible=False)
            snp_axs.set_ylim([0, 18])
            if zoomscatter_type=="brown" and c_type=="add":
                snp_axs.text(-.2, .97, r"$\bf{A}$", fontsize=13, transform=snp_axs.transAxes)
            snp_fig.savefig(f"{output_dir}/{genodata_type}_{zoomscatter_type}_{c_type}_{v_window[0]}_windows.pdf",format="pdf", bbox_inches="tight")
            plt.close(snp_fig)
    sw_lpos = []
    sw_rpos = []
    sw_pmax = []
    sw_llmax = []
    sw_spmax = []
    sw_idxmax = []
    sw_argpmax = []
    sw_rsidmax = []
    sw_chrs = []
    sw_nums = []
    sw_snps = []
    sw_raw_snps = []
    sw_raw_nums = []
    sw_type = []
    sw_ref = []
    sw_alt = []
    sw_chridxmaxs = []
    sw_genes = []
    last_col_name = r"$\hat{s}(p_{min})$" if c_type != "full" else r"$\widehat{s_1, s_2}(p_{min})$"
    for sig_window in valid_windows:
        full_window = np.arange(sig_window[0], sig_window[-1]+1)
        lpos = cdata["all_loc_per_chrom"][min(sig_window)]
        rpos = cdata["all_loc_per_chrom"][max(sig_window)]
        sw_lpos.append(lpos)
        sw_rpos.append(rpos)
        window_p_vals = cdata["all_p"][f"{c_type}_p"][full_window]
        sw_pmax.append(max(window_p_vals))
        window_argmax = np.argmax(window_p_vals)
        sw_idxmax.append(full_window[window_argmax])
        sw_argpmax.append(cdata["all_loc_per_chrom"][full_window[window_argmax]])
        sw_rsidmax.append(cdata["all_rsid"][full_window[window_argmax]])
        sw_ref.append(cdata["all_ref_allele"][full_window[window_argmax]])
        sw_alt.append(cdata["all_alt_allele"][full_window[window_argmax]])
        sw_spmax.append(cdata["all_s"][f"{c_type}_s"][full_window[np.argmax(window_p_vals)]])
        sw_llmax.append(cdata["all_ll"][f"{c_type}_ll"][full_window[np.argmax(window_p_vals)]])
        sw_chrom = int(cdata["all_chrom"][sig_window[0]])
        sw_chrs.append(sw_chrom)
        sw_nums.append(len(sig_window))
        raw_mask = ((snp_df["snp_chr"]==cdata["all_chrom"][sig_window[0]]) & (snp_df["snp_pos"]>=lpos) & (snp_df["snp_pos"]<=rpos))
        raw_snps = snp_df["snp_idx"][raw_mask].tolist()
        sw_raw_snps.append(", ".join(str(x) for x in raw_snps))
        sw_raw_nums.append(raw_mask.sum())
        sw_snps.append(", ".join(str(x) for x in sig_window))
        sw_type.append(c_type)
        sw_chridxmaxs.append(full_window[window_argmax]-cdata["all_loc"][f"chr_{sw_chrom}_idx_offset"])
        sw_genes.append("TBD")
    if len(sw_lpos) > 0:
        if c_type == "full":
            sw_spmax = [f"({s_val[0]:.4f}, {s_val[1]:.4f})" for s_val in sw_spmax]
        sw_array = np.array([sw_type, sw_chrs, sw_lpos, sw_rpos, sw_raw_nums, sw_nums, sw_pmax, sw_spmax, sw_llmax,
                             sw_idxmax, sw_argpmax, sw_rsidmax, sw_ref, sw_alt, sw_raw_snps, sw_snps, sw_chridxmaxs, sw_genes]).T
        brown_windows = pd.DataFrame(sw_array, columns=["Sel. type", "Chr.", "Left pos.", "Right pos.", "Raw", "Post", r"$-\log_{10}p_{min}$", last_col_name, r"ll at $s(p_{min})$", "SNP index of max.", "Chr pos of max.", "Lead SNP", "Ref.", "Alt.", "Raw_SNP_list", "Post_SNP_list", "SNP. index of max (on chr).", "Gene(s)"])
        brown_windows["Genomic region (hg19)"] = brown_windows.apply(combine_pos, axis=1)
        brown_windows[r"$-\log_{10}p_{min}$"] = brown_windows[r"$-\log_{10}p_{min}$"].astype(float)
        if c_type != "full":
            brown_windows[last_col_name] = brown_windows[last_col_name].astype(float)
        brown_windows = brown_windows[["Sel. type", "Chr.", "Genomic region (hg19)", "Gene(s)", "Lead SNP", "Ref.", "Alt.",
             "Raw", "Post", r"$-\log_{10}p_{min}$", last_col_name,
             r"ll at $s(p_{min})$", "SNP index of max.", "Chr pos of max.",
             "Raw_SNP_list", "Post_SNP_list", "SNP. index of max (on chr)."]]
        brown_windows.to_latex(f"{output_dir}/{genodata_type}_{c_type}_sig_windows.tex", float_format=special_format,
                               columns=["Chr.", "Genomic region (hg19)", "Gene(s)", "Lead SNP", "Ref.", "Alt.",
                                         "Raw", "Post", r"$-\log_{10}p_{min}$", last_col_name],
                               index=False, column_format="cccccccccc")

        with open(f"{output_dir}/{genodata_type}_{c_type}_sig_windows.pkl", "wb") as file:
            pickle.dump(brown_windows, file)

    #manhattan plots
    fig, axs = plt.subplots(1,1,figsize=(6.25,2.75), layout="constrained")
    axs.plot(cdata["all_loc_all_chrom"][cdata["all_chrom"]%2==1], cdata["all_p"][f"{c_type}_p"][cdata["all_chrom"]%2==1], "o", markersize=1.5, color="#888888", rasterized=True)
    axs.plot(cdata["all_loc_all_chrom"][cdata["all_chrom"]%2==0], cdata["all_p"][f"{c_type}_p"][cdata["all_chrom"]%2==0], "o", markersize=1.5, color="#87CFFF", rasterized=True)
    axs.axhline(-np.log10(snp_df["p_bh"]), ls="--", c="r", label=r"BH thresh.", lw=.75)
    axs.set_ylabel(r"$-\log_{10}(p)$")
    axs.set_xlabel("Chromosome")
    all_chr_pos_offsets = np.array([cdata["all_loc"][f"chr_{i}_pos_offset"] for i in np.arange(len(chroms)+1)+1])
    axs.xaxis.set_tick_params(length=0)
    axs.set_xticks((all_chr_pos_offsets[:-1]+all_chr_pos_offsets[1:])/2)
    chrom_labels = [str(c_i) if c_i%2 else "" for c_i in chroms]
    axs.set_xticklabels(chrom_labels)
    axs.legend()
    axs.set_ylim([0, 18])
    axs.set_xlim([0, cdata["all_loc_all_chrom"][-1]*1.001])
    fig.savefig(f"{output_dir}/{genodata_type}_{c_type}_manhattan_rasterized.pdf", format="pdf", bbox_inches="tight", dpi=600)
    plt.close(fig)

#qq plots
if "full" in classification_types:
    fig, axs = plt.subplots(1,1,figsize=(3.1, 3.1),layout="constrained")
    axins = axs.inset_axes([.67, .11, .28, .28])
    logps = [cdata['all_p'][f'{ctype}_p'] for ctype in all_classification_types]
    labels = convert_from_abbrevs(all_classification_types, shorthet=True)
    plot_qq(axs, axins, logps, labels, legend_loc="upper right", thin=True, rasterized=True)
    fig.savefig(f"{output_dir}/{genodata_type}_all_qqs_rasterized.pdf", format="pdf", bbox_inches="tight", dpi=600)