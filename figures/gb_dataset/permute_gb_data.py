import numpy as np
from emsel_util import vcf_to_useful_format
import allel
from pandas import read_csv
from tqdm import tqdm


###### MODIFY

data_dir = "data"
output_dir = "output"
genodata_type = "capture_only"
MAF_filter = .05
min_sample_filter = .1

###### DO NOT MODIFY
chroms = range(1,23)
ytg = 28.1

permute = False
rng = np.random.default_rng(0)

for chrom in tqdm(chroms):
    chrom_data_path = f"{data_dir}/GB_v54.1_{genodata_type}_c{chrom}"
    chrom_save_path = chrom_data_path + "_permuted"
    test = allel.read_vcf(chrom_data_path+".vcf")
    test_times = read_csv(f"{data_dir}/GB_v54.1_{genodata_type}_inds.table", usecols=["Genetic_ID","Date_mean"], sep="\t").to_numpy()
    test_times_perm = np.copy(test_times)
    rng.shuffle(test_times_perm[:, 1])

    geno_data = test["calldata/GT"]
    final_table = vcf_to_useful_format(test, test_times_perm, years_per_gen=ytg)

    total_fd = np.sum(final_table[:, 2::3], axis=1)
    total_ns = np.sum(final_table[:, 1::3], axis=1)

    min_fd = np.minimum(total_fd, total_ns - total_fd)
    MAF_filter_mask = min_fd > total_ns * MAF_filter

    num_samples_mask = np.sum(final_table[:, 1::3] != 0, axis=1) > 1
    anc_samples_mask = np.sum(final_table[:, 1::3], axis=1) > min_sample_filter * test["samples"].shape[0]
    combo_mask = num_samples_mask & MAF_filter_mask & anc_samples_mask

    print(f"chrom {chrom}:\nMAF filtered {final_table.shape[0] - MAF_filter_mask.sum()}\n"
          f"< 2 samples {final_table.shape[0] - num_samples_mask.sum()}\n"
          f"few samples filtered {final_table.shape[0] - anc_samples_mask.sum()}\n"
          f"overall filtered {final_table.shape[0] - combo_mask.sum()}")

    np.savetxt(chrom_save_path+".csv", final_table[combo_mask, :], delimiter="\t", fmt="%d",
               header="Each row = one replicate; each set of three columns = (sampling time, total samples, derived alleles)")