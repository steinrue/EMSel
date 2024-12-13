import numpy as np
from pathlib import Path
from tqdm import tqdm
from emsel.emsel_util import params_dict_to_str
import pickle


####### MODIFY

num_gens_list = [125]
init_dists = ["real_special"]

EM_dir = "EM/real_matched"
data_dir = "data/real_matched"


###### DO NOT MODIFY

if "matched" in EM_dir:
    EM_suffix = "Ne10496_"
    output_suffix = "real_matched_"
elif "ibdne" in EM_dir:
    EM_suffix = "Ne35313_"
    output_suffix = "ibdne_"
else:
    EM_suffix = ""
    output_suffix = ""


for num_gens in num_gens_list:
    for init_dist in init_dists:
        ndict = {}
        ndict["sel_type"] = "neutral"
        ndict["num_gens"] = num_gens
        ndict["init_dist"] = init_dist

        neutral_fname = params_dict_to_str(**ndict)
        em_path = Path(f"{EM_dir}/{neutral_fname}_{EM_suffix}EM.pkl")
        means_path = Path(f"{data_dir}/{neutral_fname}_means.txt")
        with open(em_path, "rb") as file:
            hf = pickle.load(file)
        init_mean = hf["neutral_ic"][0, :]/(hf["neutral_ic"][0, :]+hf["neutral_ic"][1, :])
        np.savetxt(means_path, np.concatenate(([hf['maf_thresh']], init_mean)), delimiter="\t", fmt="%.4f")
