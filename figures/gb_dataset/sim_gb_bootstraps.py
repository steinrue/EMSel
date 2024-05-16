import numpy as np
import pickle
from pathlib import Path
from emsel.emsel_util import get_sel_coeffs_from_type, generate_data
import matplotlib.pyplot as plt

###### MODIFY

data_dir = "data/bootstrap"
output_dir = "output"
genodata_type = "capture_only"

###### DO NOT MODIFY

pd = {"Ne": 10000,}
sampling_matrix = np.loadtxt(Path(f"{data_dir.split('/')[0]}/GB_v54.1_{genodata_type}_sample_sizes.table"), skiprows=1,dtype=int)
means_path = f"{output_dir}/GB_v54.1_{genodata_type}_means.txt"
missingness_path = f"{output_dir}/GB_v54.1_{genodata_type}_missingness.txt"

means_file = np.loadtxt(means_path, delimiter="\t")
if means_file[0] != .05:
    print("WARNING: first value of means file is not 0.05 - may not be a MAF filter")
missingness_file = np.loadtxt(missingness_path, delimiter="\t")
if missingness_file[0] != .1:
    print("WARING: first value of missingness file is not 0.1 - may not be a missingness filter")
sampling_matrix[:, 0] = sampling_matrix[-1, 0] - sampling_matrix[:, 0]
sampling_matrix = np.flip(sampling_matrix, axis=0)
num_gens = sampling_matrix[-1, 0] + 1

with open(Path(f"{output_dir.split('/')[0]}/GB_v54.1_{genodata_type}_agg_data.pkl"), "rb") as file:
    real_data_file = pickle.load(file)

temp_seed = 5

for sel_type in ["add", "dom", "het", "rec"]:
    temp_seed += 1
    with open(Path(f"{output_dir}/{genodata_type}_{sel_type}_sig_windows.pkl"), "rb") as file:
        lead_snps_table = pickle.load(file)
    pdict = {}
    pdict["sel_type"] = sel_type
    pdict["num_gens"] = num_gens
    for _, lead_snp in lead_snps_table.iterrows():
        chr = int(lead_snp["Chr."])
        chr_idx = int(lead_snp["SNP. index of max (on chr)."])
        real_data_path = Path(f"{data_dir.split('/')[0]}/GB_v54.1_{genodata_type}_c{chr}.csv")
        df = np.loadtxt(real_data_path, delimiter="\t")
        final_data = df[chr_idx, 2::3].astype(int)
        num_samples = df[chr_idx, 1::3].astype(int)
        sample_times = df[chr_idx, ::3].astype(int)
        sampling_matrix = np.vstack((sample_times, num_samples)).T
        assert sampling_matrix[:, 1].sum() == num_samples.sum()
        sel_str = float(lead_snp[r"$\hat{s}(p_{min})"])
        if sel_type == "het":
            if sel_str > 0:
                true_sel_type = "over"
            else:
                true_sel_type = "under"
                sel_str = -sel_str
        else:
            true_sel_type = sel_type
        p_init = real_data_file["all_means"][int(lead_snp["SNP index of max."])]
        rsid = lead_snp["Lead SNP"]
        pdict["sel_str"] = sel_str
        pdict["init_dist"] = p_init
        params_filename = Path(f"{data_dir}/{rsid}_{sel_type}_bootstrap_pd.pkl")
        data_filename = Path(f"{data_dir}/{rsid}_{sel_type}_data.csv")

        s1, s2 = get_sel_coeffs_from_type(true_sel_type, sel_str)
        pd["s1_true"] = s1
        pd["s2_true"] = s2
        pd["num_gens"] = num_gens
        pd["p_init"] = pdict["init_dist"]
        pd["seed"] = temp_seed
        pd["survive_only"] = True
        pd["sel_type"] = true_sel_type
        pd["small_s"] = True
        pd["normal_timestep"] = 0
        pd["num_sims"] = 1000
        pd["init_cond"] = "delta"
        pd["means_array"] = means_file
        pd["missingness_array"] = missingness_file
        pd["sampling_matrix"] = sampling_matrix
        data_dict = generate_data(pd)

        data_dict["obs_counts"] = data_dict["obs_counts"].astype(int)
        data_dict["nt"] = data_dict["nt"].astype(int)

        if data_dict["true_data"].shape[0] < pd["num_sims"]:
            continue

        pd["true_data"] = data_dict["true_data"]
        pd["p_inits"] = data_dict["p_inits"]
        fig, axs = plt.subplots(1, 1, figsize=(20, 20))
        axs.plot(data_dict["true_data"].T)
        axs.plot(np.mean(data_dict["true_data"], axis=0), color="k", lw=2)
        fig.savefig(Path(f"{output_dir}/{rsid}_{sel_type}_plot.png"), bbox_inches="tight")
        plt.close(fig)

        data_csv = np.zeros((data_dict["true_data"].shape[0], data_dict["sample_locs"].shape[0]*3))

        data_csv[:, ::3] = data_dict["sample_locs"]
        data_csv[:, 1::3] = data_dict["nt"]
        data_csv[:, 2::3] = data_dict["obs_counts"]
        with open(params_filename, "wb") as file:
            pickle.dump(pd, file)
        np.savetxt(data_filename, data_csv, delimiter="\t", fmt="%d",
                   header="Each row = one replicate; each set of three columns = (sampling time, total samples, derived alleles)")




