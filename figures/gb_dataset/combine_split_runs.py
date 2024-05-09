import pickle
from pathlib import Path
from tqdm import tqdm

###### MODIFY

EM_dir = "EM"
genodata_type = "capture_only"

###### DO NOT MODIFY

new_ugs = [["neutral", "add"], ["dom"], ["het"], ["rec"], ["full"]]
total_ugs = len(new_ugs)

chroms = range(1,23)

for chrom in tqdm(chroms):
    hmm_basepath = f"{EM_dir}/GB_v54.1_{genodata_type}_c{chrom}_EM"
    print(hmm_basepath)

    if Path(hmm_basepath + ".pkl").is_file():
        print(f"already exists: {hmm_basepath}")
        continue
    hmm_1 = Path(hmm_basepath + f"_1o{total_ugs}.pkl")

    with open(hmm_1, "rb") as file:
        hf = pickle.load(file)

    for ug_i, update_group in enumerate(new_ugs[1:], start=2):
        hmm_i = Path(hmm_basepath + f"_{ug_i}o{total_ugs}.pkl")

        with open(hmm_i, "rb") as file:
            hf_i = pickle.load(file)

        for update_type in update_group:
            hf[f"{update_type}_run"] = hf_i[f"{update_type}_run"]

        Path.unlink(hmm_i)

    with open(Path(hmm_basepath + ".pkl"), "wb") as file:
        pickle.dump(hf, file)

    Path.unlink(hmm_1)