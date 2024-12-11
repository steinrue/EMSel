import numpy as np
from pathlib import Path
import pickle
from scipy.stats import gengamma


###### MODIFY
EM_dir = "EM"
output_dir = "output"

###### DO NOT MODIFY

neutral_filename = "neutral_g125_dal_special_Ne9715"
perm_path = Path(f"{EM_dir}/{neutral_filename}_EM.pkl")
with open(perm_path, "rb") as file:
    nf = pickle.load(file)

gengamma_save_path = Path(f"{output_dir}/{neutral_filename}_gengamma_fit.pkl")


statistic = 2*(nf["full_run"]["ll_final"]-nf["neutral_ll"])

statistic[statistic <= 0] = 1e-12

lr_shift = np.median(statistic)
gg_fit = gengamma(*gengamma.fit(statistic[statistic > lr_shift] - lr_shift, floc=0, fscale=1, method="mm"))

gengamma_dict = {"lr_shift": lr_shift, "gengamma_fit": gg_fit}
with open(gengamma_save_path, "wb") as file:
    pickle.dump(gengamma_dict, file)