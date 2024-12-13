import numpy as np
from pathlib import Path
import pickle
from scipy.stats import gengamma
from emsel.emsel_util import params_dict_to_str

###### MODIFY
num_gens_list = [101, 251, 1001]
init_dists = [.005, .25, "recip"]

#num_gens_list = [125]
#init_dists = ["real_special"]

EM_dir = "EM/pure_sim"
output_dir = "output/pure_sim"
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

        neutral_filename = params_dict_to_str(**ndict)

        resim_path = Path(f"{EM_dir}/{neutral_filename}_resim_{EM_suffix}EM.pkl")
        with open(resim_path, "rb") as file:
            nf = pickle.load(file)

        resim_gengamma_save_path = Path(f"{output_dir}/{neutral_filename}_{output_suffix}gengamma_resim_fit.pkl")


        statistic = 2*(nf["full_run"]["ll_final"]-nf["neutral_ll"])
        #print(np.min(statistic))

        statistic[statistic <= 0] = 1e-12

        lr_shift = np.median(statistic)
        gg_fit = gengamma(*gengamma.fit(statistic[statistic > lr_shift] - lr_shift, floc=0, fscale=1, method="mm"))

        gengamma_dict = {"lr_shift": lr_shift, "gengamma_fit": gg_fit}
        with open(resim_gengamma_save_path, "wb") as file:
            pickle.dump(gengamma_dict, file)