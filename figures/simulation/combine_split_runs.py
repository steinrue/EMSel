import pickle
from pathlib import Path
from tqdm import tqdm

###### MODIFY

EM_dir = "EM"

###### DO NOT MODIFY

new_ugs = [["neutral", "add"], ["dom"], ["het"], ["rec"], ["full"]]
total_ugs = len(new_ugs)

num_files = 0
for fpath in Path(EM_dir).iterdir():
    if fpath.name.rpartition("_")[-1] == f"1o{total_ugs}.pkl":
        hmm_basepath = str(fpath).rpartition("_")[0]
        print(hmm_basepath)

        if Path(hmm_basepath + ".pkl").is_file():
            print(f"already exists: {hmm_basepath}")
            continue

        num_files += 1

        with open(fpath, "rb") as file:
            hf = pickle.load(file)

        for ug_i, update_group in enumerate(new_ugs[1:], start=2):
            hmm_i = Path(hmm_basepath + f"_{ug_i}o{total_ugs}.pkl")

            if hmm_i.is_file():
                with open(hmm_i, "rb") as file:
                    hf_i = pickle.load(file)
                for update_type in update_group:
                    hf[f"{update_type}_run"] = hf_i[f"{update_type}_run"]
                Path.unlink(hmm_i)
            else:
                print(f"No file to combine - {hmm_i.name}")
        with open(Path(hmm_basepath + ".pkl"), "wb") as file:
            pickle.dump(hf, file)

        Path.unlink(fpath)
print(f"{num_files} files combined!")