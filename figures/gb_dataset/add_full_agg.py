import numpy as np
from pathlib import Path
import pickle
from emsel.emsel_util import bh_correct, correct_for_het, get_1d_s_data_from_type, get_llg_array, get_llgka_array, bh_procedure_2, classify_full_run
from copy import deepcopy
from scipy.stats import chi2

###### MODIFY

data_dir = "data"
EM_dir = "EM"
output_dir = "output"
genodata_type = "capture_only"

###### DO NOT MODIFY
def get_maxbased_llr(hf, onep_types, uncon=False):
    temp_lls = np.zeros_like(hf["neutral_ll"]) - np.inf
    for onep_run in onep_types:
        temp_lls = np.maximum(temp_lls, hf[f"{onep_run}_run"]["ll_final"])
    if uncon:
        temp_lls = np.maximum(temp_lls, hf["full_run"]["ll_final"]-1)
    return 2*(temp_lls-hf["neutral_ll"])

chroms = range(1,23)
alpha = .05

onep_types = ["add", "dom", "rec", "het"]
full_p = np.array([0])
all_llg_array = np.zeros((1,len(onep_types) + 4))
all_llr_array = []
all_het_s_data = []

complete_agg_data_path = Path(f"{output_dir}/GB_v54.1_{genodata_type}_agg_data.pkl")
for chrom in chroms:
    base_data_path = Path(f"{data_dir}/GB_v54.1_{genodata_type}_c{chrom}")

    #aggregating various things
    em_path = Path(f"{EM_dir}/GB_v54.1_{genodata_type}_c{chrom}_EM.pkl")
    with open(em_path, "rb") as file:
        hf = pickle.load(file)


    temp_llg_array = get_llg_array(hf, onep_types, classify_full_run(hf["full_run"]["s_final"])[0])
    all_llg_array = np.vstack((all_llg_array, temp_llg_array))
    temp_hs_data = hf["het_run"]["s_final"][0, :][:, np.newaxis]
    all_het_s_data.append(temp_hs_data)
    temp_llr_data = get_maxbased_llr(hf, onep_types, uncon=False)
    all_llr_array.append(temp_llr_data)

all_het_s_data = np.concatenate(all_het_s_data)
all_llr_array = np.concatenate(all_llr_array)

with open(complete_agg_data_path, "rb") as file:
    cdata = pickle.load(file)

all_llg_array = all_llg_array[1:, :]
all_llgka_array = get_llgka_array(all_llg_array[:, :, np.newaxis], k=np.inf, alpha=0)
p_full_bh, full_p_vals, full_logp_vals, full_classes = bh_procedure_2([all_llr_array], [all_llgka_array], chi2(1), 0, alpha, bh=True)


final_classified_vals = correct_for_het(deepcopy(full_classes), all_het_s_data, onep_types)
full_p_vals = full_logp_vals[0]
full_classes = final_classified_vals[0].flatten()
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

with open(complete_agg_data_path, "wb") as file:
    pickle.dump(cdata, file)