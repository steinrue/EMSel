from pathlib import Path

###### MODIFY
EM_dir = Path('EM/bootstrap')
data_dir = Path('data/bootstrap')
genodata_type = "capture_only"

qsub_dir = Path('qsubs')
max_run_hours = 20
num_cores = 30
mem_required = 30

Ne = 9715


#this directory must be made beforehand


###### DO NOT MODIFY

chroms = range(1, 23)
new_ugs = [["neutral", "add"], ["dom"], ["het"], ["rec"], ["full"]]
total_ugs = len(new_ugs)

def writeQsubs():
    hmm_cmds = []

    # get a meta file going
    script_file = Path("meta_gb_EM.sh")
    with open(script_file, "w") as file:
        for fpath in data_dir.iterdir():
            if "bootstrap" not in data_dir.name and genodata_type not in fpath.name:
                continue
            for ug_i, update_group in enumerate(new_ugs):
                if fpath.suffix == ".vcf" or fpath.suffix == ".csv":
                    out_name = str(Path(EM_dir, f"{fpath.stem}_EM_{ug_i+1}o{total_ugs}"))
                    sbatchfile = Path(qsub_dir, f"{fpath.stem}_EM_{ug_i+1}o{total_ugs}.sbatch")
                    sbatchOutFile = Path(qsub_dir, f"{fpath.stem}_EM_{ug_i+1}o{total_ugs}.sout")
                    sbatchErrFile = Path(qsub_dir, f"{fpath.stem}_EM_{ug_i+1}o{total_ugs}.serr")
                    # put together a nice command?
                    if fpath.suffix == ".vcf" and "only" in genodata_type:
                        hmm_cmd = f"emsel {fpath} {out_name} --time_before_present {'--save_csv' if ug_i == 0 else ''} --info_file {data_dir}/GB_v54.1_{genodata_type}_inds.table --info_cols Genetic_ID Date_mean -ytg 28.1 --full_output --num_cores {num_cores} --selection_modes {' '.join(u_i for u_i in update_group)} --no_neutral --progressbar -Ne {Ne}"
                    else:
                        if ug_i > 0:
                            continue
                        vcf_ver = Path(fpath).with_suffix(".vcf")
                        if fpath.suffix == ".vcf" and "SG" in genodata_type:
                            hmm_cmd = f"emsel {fpath} {out_name} --time_before_present --save_csv --info_file {data_dir}/GB_v54.1_{genodata_type}_inds.table --info_cols Genetic_ID Date_mean -ytg 28.1 --full_output --num_cores {num_cores} --selection_modes neutral add --progressbar -Ne {Ne}"
                        elif vcf_ver.is_file():
                            continue
                        if "bootstrap" in data_dir.name:
                            sel_type = fpath.stem.split("_")[1]
                            hmm_cmd = f"emsel {fpath} {out_name} --num_cores {num_cores} --time_after_zero --full_output --selection_modes {sel_type} --no_neutral --progressbar -Ne {Ne}"
                        elif "permuted" in fpath.name:
                            if "100k" in fpath.name:
                                hmm_cmd = f"emsel {fpath} {out_name} --num_cores {num_cores} --time_after_zero --full_output --selection_modes neutral add full --progressbar -Ne {Ne}"
                            else:
                                hmm_cmd = f"emsel {fpath} {out_name} --num_cores {num_cores} --time_after_zero --full_output --selection_modes neutral add --progressbar -Ne {Ne}"
                    if Path(out_name+'.csv').is_file() and Path(out_name+".csv").stat().st_size > 0:
                        print(f"File already exists: {out_name}! continuing.")
                        continue

                    all_path = Path(out_name.rpartition("_")[0] + ".pkl")
                    if (all_path.is_file() and all_path.stat().st_size > 0):
                        print(f"File already exists: {all_path.name}! continuing.")
                        continue
                    if hmm_cmd not in hmm_cmds:
                        # write a script for the cluster
                        with open(sbatchfile, 'w') as qsubScriptFile:
                            qsubScriptFile.write("#!/bin/bash\n" + \
                                                 f"#SBATCH -J {fpath.name.rpartition('_')[-1]}\n" + \
                                                 f"#SBATCH --time={max_run_hours}:00:00\n" + \
                                                 f"#SBATCH --cpus-per-task={num_cores}\n" + \
                                                 f"#SBATCH --mem={mem_required}gb\n" + \
                                                 f"#SBATCH -o {sbatchOutFile}\n" + \
                                                 f"#SBATCH -e {sbatchErrFile}\n" + \
                                                 f"{hmm_cmd}\n")
                        # and add to meta file
                        file.write(f"sbatch {sbatchfile}\n")
                        print(hmm_cmd)
                        hmm_cmds.append(hmm_cmd)
    print(f"{len(hmm_cmds)} files to be run!")

def main():
    writeQsubs()


if __name__ == "__main__":
    main()