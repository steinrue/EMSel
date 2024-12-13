import numpy as np
from pathlib import Path
import pickle
from scipy.stats import chi2
from emsel.emsel_util import bh_correct, get_1d_s_data_from_type

###### MODIFY

data_dir = "data"
EM_dir = "EM"
output_dir = "output"
genodata_type = "capture_only"
suffix = ""
classification_types = ["add", "dom", "het", "rec", "full"]

###### DO NOT MODIFY
chroms = range(1,23)

alpha = .05

all_s = {}
all_p = {}
all_ll = {}
for c_type in classification_types:
    all_p[f"{c_type}_p"] = np.array([0])
    all_s[f"{c_type}_s"] = np.array([0]) if c_type != "full" else np.zeros((1,2))
    all_ll[f"{c_type}_ll"] = np.array([0])
all_ll["neutral_ll"] = np.array([0])
all_loc = {}
all_maf = {}
all_loc["chr_1_idx_offset"] = 0
all_loc["chr_1_pos_offset"] = 0
all_rsid = np.array([0])
all_chrom = np.array([0])
all_loc_all_chrom = np.array([0], dtype=np.int64)
all_loc_per_chrom = np.array([0])
all_fd = np.zeros((1, 5))
all_ns = np.zeros((1, 5))
all_st = np.zeros((1, 5))
all_ref_allele = np.array([0])
all_alt_allele = np.array([0])
all_missingness = np.zeros(1)
all_means = np.zeros(1)
start_pos = 0

num_bins = 5
superfinal_binned_data = np.zeros((1, num_bins * 3))

binned_path = Path(f"{data_dir}/GB_v54.1_{genodata_type}_complete_data_binned.csv")
complete_agg_data_path = Path(f"{output_dir}/GB_v54.1_{genodata_type}_{suffix}agg_data.pkl")
means_path = Path(f"{output_dir}/GB_v54.1_{genodata_type}_means.txt")
missingness_path = Path(f"{output_dir}/GB_v54.1_{genodata_type}_missingness.txt")
for chrom in chroms:
    agg_data = {}
    base_data_path = Path(f"{data_dir}/GB_v54.1_{genodata_type}_c{chrom}")
    df = np.loadtxt(base_data_path.with_suffix(base_data_path.suffix+".csv"), delimiter="\t")
    final_data = df[:, 2::3].astype(int)
    num_samples = df[:, 1::3].astype(int)
    sample_times = df[:, ::3].astype(int)

    #aggregating various things
    em_path = Path(f"{EM_dir}/GB_v54.1_{genodata_type}_c{chrom}_{suffix}EM.pkl")
    with open(em_path, "rb") as file:
        hf = pickle.load(file)
    init_mean = hf["neutral_ic"][0, :]/(hf["neutral_ic"][0, :]+hf["neutral_ic"][1, :])
    agg_data["neut_init_mean"] = init_mean
    neutral_ll = hf["neutral_ll"]
    pos = hf["pos"]
    a_ns = np.sum(num_samples, axis=1)
    a_fd = np.sum(final_data, axis=1)
    agg_data["pos"] = pos
    agg_data["snp_ids"] = hf["snp_ids"]
    agg_data["a_numsamples"] = a_ns
    agg_data["a_freq"] = a_fd / a_ns
    agg_data["a_miss"] = 1-(a_ns / hf["max_samples"])
    agg_data["ref_allele"] = hf["ref_allele"]
    agg_data["alt_allele"] = hf["alt_allele"]
    agg_data["neutral_ll"] = neutral_ll
    for classification_type in classification_types:
        if classification_type == "full":
            # to be overwritten by add_full_agg.py
            agg_data[f"{classification_type}_ll_vals"] = np.zeros_like(agg_data["pos"])
            agg_data[f"{classification_type}_p_vals"] = np.zeros_like(agg_data["pos"])
            agg_data[f"{classification_type}_s_vals"] = hf[f"{classification_type}_run"]["s_final"].T
        else:
            run_ll = hf[f"{classification_type}_run"]["ll_final"]
            llr_real = 2 * (run_ll - neutral_ll)
            llr = np.copy(llr_real)
            llr[llr <= 0] = 1e-12
            p_vals = -chi2(1).logsf(llr)/np.log(10)
            agg_data[f"{classification_type}_ll_vals"] = llr
            agg_data[f"{classification_type}_p_vals"] = p_vals
            agg_data[f"{classification_type}_itercount"] = hf[f"{classification_type}_run"]["itercount_hist"]
            agg_data[f"{classification_type}_s_vals"] = get_1d_s_data_from_type(hf[f"{classification_type}_run"]["s_final"], classification_type)
            if not np.all(np.isfinite(p_vals)):
                print(f"c {chrom} {classification_type} illegal stuff happening")
    for c_type in classification_types:
        all_p[f"{c_type}_p"] = np.concatenate((all_p[f"{c_type}_p"], agg_data[f"{c_type}_p_vals"]))
        all_s[f"{c_type}_s"] = np.concatenate(
            (all_s[f"{c_type}_s"], agg_data[f"{c_type}_s_vals"]))
        all_ll[f"{c_type}_ll"] = np.concatenate(
            (all_ll[f"{c_type}_ll"], agg_data[f"{c_type}_ll_vals"]))
        if c_type == classification_types[0]:
            all_loc_all_chrom = np.concatenate((all_loc_all_chrom, agg_data["pos"].astype(np.int64)+start_pos))
            all_loc_per_chrom = np.concatenate((all_loc_per_chrom, agg_data["pos"]))
            start_pos += np.max(agg_data["pos"])
            all_loc[f"chr_{chrom}"] = agg_data["pos"]
            all_loc[f"chr_{chrom+1}_idx_offset"] = all_loc[f"chr_{chrom}_idx_offset"] + all_loc[f"chr_{chrom}"].shape[0]
            all_loc[f"chr_{chrom+1}_pos_offset"] = all_loc[f"chr_{chrom}_pos_offset"] + all_loc[f"chr_{chrom}"][-1]
            all_missingness = np.concatenate((all_missingness, agg_data["a_miss"]))
            all_means = np.concatenate((all_means, agg_data["neut_init_mean"]))
            all_rsid = np.concatenate((all_rsid, agg_data["snp_ids"]))
            all_chrom = np.concatenate((all_chrom, np.zeros_like(agg_data["pos"])+chrom))
            all_ref_allele = np.concatenate((all_ref_allele, agg_data["ref_allele"]))
            all_alt_allele = np.concatenate((all_alt_allele, agg_data["alt_allele"]))
            all_ll["neutral_ll"] = np.concatenate((all_ll["neutral_ll"], agg_data["neutral_ll"]))
            if chrom == 1:
                all_rsid = all_rsid[1:]
                all_chrom = all_chrom[1:]
                all_loc_all_chrom = all_loc_all_chrom[1:]
                all_loc_per_chrom = all_loc_per_chrom[1:]
                all_missingness = all_missingness[1:]
                all_means = all_means[1:]
                all_ref_allele = all_ref_allele[1:]
                all_alt_allele = all_alt_allele[1:]
                all_ll["neutral_ll"] = all_ll["neutral_ll"][1:]
                MAF_THRESH = hf["maf_thresh"]
                MISSING_THRESH = hf["min_sample_density"]
        if chrom == 1:
            all_p[f"{c_type}_p"] = all_p[f"{c_type}_p"][1:]
            all_s[f"{c_type}_s"] = all_s[f"{c_type}_s"][1:]
            all_ll[f"{c_type}_ll"] = all_ll[f"{c_type}_ll"][1:]


    #binning the data for trajectories
    if chrom == 1:
        min_c1_st = np.min(sample_times)
        max_c1_st = np.max(sample_times)
    else:
        assert min_c1_st == np.min(sample_times)
        assert max_c1_st == np.max(sample_times)
    cutoffs = np.linspace(min_c1_st - 1, max_c1_st + 1, endpoint=True, num=num_bins + 1)
    idxs = np.searchsorted(sample_times[0, :], cutoffs)
    binned_num_samples = np.zeros((num_samples.shape[0], cutoffs.shape[0] - 1))
    binned_final_data = np.zeros((num_samples.shape[0], cutoffs.shape[0] - 1))
    binned_sample_times = np.zeros((num_samples.shape[0], cutoffs.shape[0] - 1))
    for i in range(cutoffs.shape[0] - 1):
        binned_num_samples[:, i] = np.sum(num_samples[:, idxs[i]:idxs[i + 1]], axis=1)
        binned_final_data[:, i] = np.sum(final_data[:, idxs[i]:idxs[i + 1]], axis=1)
        binned_sample_times[:, i] = (cutoffs[i] + cutoffs[i + 1]) / 2

    binned_csv = np.zeros((binned_final_data.shape[0], binned_final_data.shape[1] * 3))
    binned_csv[:, ::3] = binned_sample_times
    binned_csv[:, 1::3] = binned_num_samples
    binned_csv[:, 2::3] = binned_final_data
    superfinal_binned_data = np.vstack((superfinal_binned_data, binned_csv))

for classification_type in classification_types:
    if classification_type == "full":
        continue
    data_to_save = {}
    p_bh, bh_locs = bh_correct(np.power(10, -all_p[f"{classification_type}_p"]), alpha, yekutieli=False)
    data_to_save["snp_idx"] = bh_locs
    data_to_save["snp_pos"] = all_loc_per_chrom[bh_locs]
    data_to_save["snp_rsid"] = all_rsid[bh_locs]
    data_to_save["snp_chr"] = all_chrom[bh_locs]
    data_to_save["p"] = all_p[f"{classification_type}_p"][bh_locs]
    data_to_save["p_bh"] = p_bh

    with open(Path(f"{output_dir}/GB_v54.1_{genodata_type}_{classification_type}_bh.pkl"), "wb") as file:
        pickle.dump(data_to_save, file)

cdata = {}
cdata["all_p"] = all_p
cdata["all_s"] = all_s
cdata["all_ll"] = all_ll
cdata["all_loc_all_chrom"] = all_loc_all_chrom
cdata["all_loc_per_chrom"] = all_loc_per_chrom
cdata["all_chrom"] = all_chrom
cdata["all_rsid"] = all_rsid
cdata["all_loc"] = all_loc
cdata["all_missingness"] = all_missingness
cdata["all_means"] = all_means
cdata["all_ref_allele"] = all_ref_allele
cdata["all_alt_allele"] = all_alt_allele
cdata["maf_thresh"] = MAF_THRESH
cdata["missing_thresh"] = MISSING_THRESH

with open(complete_agg_data_path, "wb") as file:
    pickle.dump(cdata, file)

superfinal_binned_data = superfinal_binned_data[1:, :]
np.savetxt(binned_path, superfinal_binned_data, delimiter="\t", fmt="%d",
           header="Each row = one replicate; each set of three columns = (sampling time, total samples, derived alleles)")

np.savetxt(means_path, np.concatenate(([MAF_THRESH], cdata["all_means"])), delimiter="\t", fmt="%.4f")
np.savetxt(missingness_path, np.concatenate(([MISSING_THRESH], cdata["all_missingness"])), delimiter="\t", fmt="%.4f")