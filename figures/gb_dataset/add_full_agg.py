import numpy as np
from pathlib import Path
import pickle
from emsel.emsel_util import bh_correct, get_1d_s_data_from_type, get_llg_array, get_llgka_array, full_bh_procedure, classify_full_run

###### MODIFY

data_dir = "data"
EM_dir = "EM"
output_dir = "output"
genodata_type = "capture_only"

###### DO NOT MODIFY
chroms = range(1,23)
alpha = .05

onep_types = ["add", "dom", "rec"]
full_p = np.array([0])
all_llg_array = np.zeros((1,len(onep_types) + 4))

binned_path = Path(f"{data_dir}/GB_v54.1_{genodata_type}_complete_data_binned.csv")
complete_agg_data_path = Path(f"{output_dir}/GB_v54.1_{genodata_type}_agg_data.pkl")
means_path = Path(f"{output_dir}/GB_v54.1_{genodata_type}_means.txt")
missingness_path = Path(f"{output_dir}/GB_v54.1_{genodata_type}_missingness.txt")
gengamma_path = Path(f"{data_dir}/gengamma_params.pkl")
for chrom in chroms:
    base_data_path = Path(f"{data_dir}/GB_v54.1_{genodata_type}_c{chrom}")

    #aggregating various things
    em_path = Path(f"{EM_dir}/GB_v54.1_{genodata_type}_c{chrom}_EM.pkl")
    with open(em_path, "rb") as file:
        hf = pickle.load(file)
    temp_llg_array = get_llg_array(hf, onep_types, classify_full_run(hf["full_run"]["s_final"])[0])
    all_llg_array = np.vstack((all_llg_array, temp_llg_array))

with open(complete_agg_data_path, "rb") as file:
    cdata = pickle.load(file)

with open(gengamma_path, "rb") as file:
    gengamma_dict = pickle.load(file)

all_llg_array = all_llg_array[1:, :]
all_llgka_array = get_llgka_array(all_llg_array[:, :, np.newaxis], k=gengamma_dict["k_opt"], alpha=0)
p_full_bh, full_p_vals, full_classes = full_bh_procedure([all_llgka_array], gengamma_dict["gengamma_fit"],
                                                         gengamma_dict["lr_shift"], alpha, bh=True)
full_p_vals = -np.log10(full_p_vals[0].flatten())
full_classes = full_classes[0].flatten()
full_data_to_save = {}
bh_locs = np.where(full_classes > 0)[0]
full_data_to_save["snp_idx"] = bh_locs
full_data_to_save["snp_pos"] = cdata["all_loc_per_chrom"][bh_locs]
full_data_to_save["snp_rsid"] = cdata["all_rsid"][bh_locs]
full_data_to_save["snp_chr"] = cdata["all_chrom"][bh_locs]
full_data_to_save["p"] = full_p_vals[bh_locs]
full_data_to_save["p_bh"] = p_full_bh
full_data_to_save["classes"] = full_classes[bh_locs]

with open(Path(f"{output_dir}/GB_v54.1_{genodata_type}_full_bh.pkl"),"wb") as file:
    pickle.dump(full_data_to_save, file)


cdata["all_p"]["full_p"] = full_p_vals

#this is only necessary so that a later script doesn't throw an error.
cdata["all_s"]["full_s"] = cdata["all_s"]["add_s"]

with open(complete_agg_data_path, "wb") as file:
    pickle.dump(cdata, file)
