import numpy as np
from pathlib import Path
from tqdm import tqdm
from emsel.emsel_util import params_dict_to_str
import pickle


####### MODIFY

EM_dir = "EM"
data_dir = "data"

###### DO NOT MODIFY
init_mean = np.array([0])
suffixes = ["uncon", "neut"]
em_path = Path(f"{EM_dir}/horse_data_full_EM.pkl")
with open(em_path, "rb") as file:
    hf = pickle.load(file)
for r_i, run in enumerate([hf["full_run"]["ic_dist"], hf["neutral_ic"]]):
    init_mean_temp = run[0, :]/(run[0, :]+run[1, :])
    init_mean_2 = np.concatenate((init_mean, [init_mean_temp.flatten()[0]]))
    np.savetxt(f"{data_dir}/full_{suffixes[r_i]}_means.txt", init_mean_2, delimiter="\t", fmt="%.4f")