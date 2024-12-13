import numpy as np
from pathlib import Path
import pickle
from scipy.stats import chi2
from scipy.interpolate import CubicSpline
from emsel.emsel_util import bh_correct, get_1d_s_data_from_type
from tqdm import tqdm
import matplotlib.pyplot as plt

###### MODIFY
EM_dir = "EM/real_matched"
###### DO NOT MODIFY

Nes = np.geomspace(500,500000, 21, dtype=int)

condo_lls = np.zeros_like(Nes, dtype=float)
for Ne_i, Ne in enumerate(Nes):
        unif_em_path = Path(f"{EM_dir}/neutral_g125_dal_special_9734plzz_HMMNe{Ne}_cond_unif_EM.pkl")
        with open(unif_em_path, "rb") as file:
            unif_hf = pickle.load(file)

        uncondo_sum = unif_hf["neutral_ll"].sum()
        condo_sum = uncondo_sum - unif_hf["cond_correction_ll"].sum()
        condo_lls[Ne_i] += condo_sum


Nes_space = np.geomspace(500, 500000, 5000)
condo_spline = CubicSpline(Nes, condo_lls)
condo_spline_output = condo_spline(Nes_space)

print(f"{EM_dir.split('/')[-1]}: Estimated Ne: {Nes_space[np.argmax(condo_spline_output)]:.0f}")