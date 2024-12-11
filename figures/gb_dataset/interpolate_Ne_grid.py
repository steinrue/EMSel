import numpy as np
from pathlib import Path
import pickle
from scipy.stats import chi2
from scipy.interpolate import CubicSpline
from emsel.emsel_util import bh_correct, get_1d_s_data_from_type
from tqdm import tqdm
import matplotlib.pyplot as plt

###### MODIFY

data_dir = "data"
EM_dir = "EM/Ne"
output_dir = "output"
genodata_type = "capture_only"

###### DO NOT MODIFY

Nes = np.geomspace(500,500000, 21, dtype=int)
chroms = range(1,23)

midfix = "cond"
condo_lls = np.zeros_like(Nes, dtype=float)
for Ne_i, Ne in enumerate(Nes):
    for chrom in chroms:
        unif_em_path = Path(f"{EM_dir}/GB_v54.1_{genodata_type}_c{chrom}_{midfix}_Ne{Ne}_unif_EM.pkl")
        with open(unif_em_path, "rb") as file:
            unif_hf = pickle.load(file)

        uncondo_sum = unif_hf["neutral_ll"].sum()
        condo_sum = uncondo_sum - unif_hf["cond_correction_ll"].sum()
        condo_lls[Ne_i] += condo_sum


Nes_space = np.geomspace(500, 500000, 5000)
condo_spline = CubicSpline(Nes, condo_lls)
condo_spline_output = condo_spline(Nes_space)

print(f"Estimated Ne: {Nes_space[np.argmax(condo_spline_output)]:.0f}")