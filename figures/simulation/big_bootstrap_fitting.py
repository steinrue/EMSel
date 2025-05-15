import numpy as np
from pathlib import Path
import pickle
from scipy.stats import gengamma, chi2
from emsel.emsel_util import params_dict_to_str, get_llg_array, get_llgka_array

###### MODIFY
num_gens_list = [101, 251, 1001]
init_dists = [.005, .25, "recip"]

num_gens_list = [125]
init_dists = ["real_special"]

EM_dir = "EM/ibdne/big"
output_dir = "output/ibdne"
onep_types = ["add", "dom", "rec", "het"]
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

        n_path = Path(f"{EM_dir}/{neutral_filename}_fit_{EM_suffix}EM.pkl")
        with open(n_path, "rb") as file:
            nf = pickle.load(file)


        #resim_gengamma_save_path = Path(f"{output_dir}/{neutral_filename}_{output_suffix}gengamma_resim_fit.pkl")
        #resim_chi2_save_path = Path(f"{output_dir}/{neutral_filename}_{output_suffix}chi2_resim_fit.pkl")
        final_gengamma_save_path = Path(f"{output_dir}/{neutral_filename}_{output_suffix}final_gengamma_fit.pkl")
        finaluncon_gengamma_save_path = Path(f"{output_dir}/{neutral_filename}_{output_suffix}finaluncon_gengamma_fit.pkl")
        final_chi2_save_path = Path(f"{output_dir}/{neutral_filename}_{output_suffix}final_chi2_fit.pkl")
        finaluncon_chi2_save_path = Path(f"{output_dir}/{neutral_filename}_{output_suffix}finaluncon_chi2_fit.pkl")

        #
        # statistic = 2*(nf["full_run"]["ll_final"]-nf["neutral_ll"])
        # #print(np.min(statistic))
        #
        # statistic[statistic <= 0] = 1e-12
        #
        # lr_shift = np.median(statistic)
        # gg_fit = gengamma(*gengamma.fit(statistic[statistic > lr_shift] - lr_shift, floc=0, fscale=1, method="mm"))
        # chi2_fit = chi2(*chi2.fit(statistic[statistic > lr_shift]-lr_shift, floc=0, fscale=1, method="mm"))
        #
        # gengamma_dict = {"lr_shift": lr_shift, "dist_fit": gg_fit}
        # with open(resim_gengamma_save_path, "wb") as file:
        #     pickle.dump(gengamma_dict, file)
        #
        # chi2_dict = {"lr_shift": lr_shift, "dist_fit": chi2_fit}
        # with open(resim_chi2_save_path, "wb") as file:
        #     pickle.dump(chi2_dict, file)
        #
        # raise Error


        neutral_llg = get_llg_array(nf, onep_types, np.zeros_like(nf["neutral_ll"]))
        neutral_llg[:, len(onep_types)+1:] -= 1

        final_statistic = 2*(np.max(neutral_llg[:, 1:len(onep_types)+1], axis=1) - nf["neutral_ll"])
        final_statistic[final_statistic <= 0] = 1e-12
        final_lr_shift = np.median(final_statistic)
        final_gg_fit = gengamma(*gengamma.fit(final_statistic[final_statistic > final_lr_shift] - final_lr_shift, floc=0, fscale=1, method="mm"))
        final_chi2_fit = chi2(*chi2.fit(final_statistic[final_statistic > final_lr_shift] - final_lr_shift, floc=0, fscale=1, method="mm"))

        final_gengamma_dict = {"lr_shift": final_lr_shift, "dist_fit": final_gg_fit}
        with open(final_gengamma_save_path, "wb") as file:
            pickle.dump(final_gengamma_dict, file)

        final_chi2_dict = {"lr_shift": final_lr_shift, "dist_fit": final_chi2_fit}
        with open(final_chi2_save_path, "wb") as file:
            pickle.dump(final_chi2_dict, file)

        finaluncon_statistic = 2*(np.max(neutral_llg[:, 1:], axis=1) - nf["neutral_ll"])
        finaluncon_statistic[finaluncon_statistic <= 0] = 1e-12
        finaluncon_lr_shift = np.median(finaluncon_statistic)
        finaluncon_gg_fit = gengamma(*gengamma.fit(finaluncon_statistic[finaluncon_statistic > finaluncon_lr_shift] - finaluncon_lr_shift, floc=0, fscale=1, method="mm"))
        finaluncon_chi2_fit = chi2(*chi2.fit(finaluncon_statistic[finaluncon_statistic > finaluncon_lr_shift] - finaluncon_lr_shift, floc=0, fscale=1, method="mm"))

        finaluncon_gengamma_dict = {"lr_shift": finaluncon_lr_shift, "dist_fit": finaluncon_gg_fit}
        with open(finaluncon_gengamma_save_path, "wb") as file:
            pickle.dump(finaluncon_gengamma_dict, file)

        finaluncon_chi2_dict = {"lr_shift": finaluncon_lr_shift, "dist_fit": finaluncon_chi2_fit}
        with open(finaluncon_chi2_save_path, "wb") as file:
            pickle.dump(finaluncon_chi2_dict, file)