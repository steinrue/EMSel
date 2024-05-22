from pathlib import Path

###### MODIFY
EM_dir = Path('EM/param_variation')
data_dir = Path('data/param_variation')

qsub_dir = Path('qsubs')
max_run_hours = 10
num_cores = 30
huge_mem_required = 50
large_mem_required = 20
small_mem_required = 8


###### DO NOT MODIFY
new_ugs = [["neutral", "add"], ["dom"], ["het"], ["rec"], ["full"]]
total_ugs = len(new_ugs)

def writeQsubs():
    hmm_cmds = []
    # get a meta file going
    script_file = Path("meta_sim_EM.sh")
    with open(script_file, "w") as file:
        for fpath in data_dir.iterdir():
            if fpath.suffix == ".csv":
                for ug_i, update_group in enumerate(new_ugs):
                    if "dal_special" in fpath.name:
                        mem = large_mem_required
                        out_name = str(Path(EM_dir, f"{fpath.name.rpartition('_')[0]}_EM_{ug_i+1}o{total_ugs}"))
                        sbatchfile = Path(qsub_dir, f"{fpath.name.rpartition('_')[0]}_EM_{ug_i+1}o{total_ugs}.sbatch")
                        sbatchOutFile = Path(qsub_dir, f"{fpath.name.rpartition('_')[0]}_EM_{ug_i+1}o{total_ugs}.sout")
                        sbatchErrFile = Path(qsub_dir, f"{fpath.name.rpartition('_')[0]}_EM_{ug_i+1}o{total_ugs}.serr")
                        if "perm" in fpath.name:
                            if ug_i > 0:
                                continue
                            hmm_cmd = f"emsel {fpath} {out_name} --time_after_zero --full_output --num_cores {num_cores} --selection_modes neutral add --progressbar"
                        else:
                            hmm_cmd = f"emsel {fpath} {out_name} --time_after_zero --full_output --num_cores {num_cores} --selection_modes {' '.join(u_i for u_i in update_group)} --no_neutral --progressbar"

                        all_path = Path(out_name.rpartition("_")[0] + ".pkl")
                        if (all_path.is_file() and all_path.stat().st_size > 0):
                            print(f"File already exists: {all_path.name}! continuing.")
                            continue
                    else:
                        mem = small_mem_required
                        #don't need to split up the non-data-matched simulations, they're smaller
                        if ug_i > 0:
                            continue
                        out_name = str(Path(EM_dir, f"{fpath.name.rpartition('_')[0]}_EM"))
                        sbatchfile = Path(qsub_dir, f"{fpath.name.rpartition('_')[0]}_EM.sbatch")
                        sbatchOutFile = Path(qsub_dir, f"{fpath.name.rpartition('_')[0]}_EM.sout")
                        sbatchErrFile = Path(qsub_dir, f"{fpath.name.rpartition('_')[0]}_EM.serr")
                        if "param_variation" in data_dir.name:
                            sel_mode = fpath.name.split("_")[0]
                            hmm_cmd = f"emsel {fpath} {out_name} --time_after_zero --full_output -maf 0 --min_sample_density 0 --num_cores {num_cores} --selection_modes {sel_mode} --no_neutral --progressbar"
                        else:
                            hmm_cmd = f"emsel {fpath} {out_name} --time_after_zero --full_output -maf 0 --min_sample_density 0 --num_cores {num_cores} --selection_modes all --progressbar"
                    if (Path(out_name).with_suffix(".csv").is_file() and Path(out_name).with_suffix(".csv").stat().st_size > 0):
                        print(f"File already exists: {out_name}! continuing.")
                        continue

                    if hmm_cmd not in hmm_cmds:
                        # write a script for the cluster
                        with open(sbatchfile, 'w') as qsubScriptFile:
                            qsubScriptFile.write("#!/bin/bash\n" + \
                                                 f"#SBATCH -J {fpath.name[:5]}\n" + \
                                                 f"#SBATCH --time={max_run_hours}:00:00\n" + \
                                                 f"#SBATCH --cpus-per-task={num_cores}\n" + \
                                                 f"#SBATCH --mem={mem}gb\n" + \
                                                 f"#SBATCH -o {sbatchOutFile}\n" + \
                                                 f"#SBATCH -e {sbatchErrFile}\n" + \
                                                 f"{hmm_cmd}\n")
                        # and add to meta file
                        file.write(f"sbatch {sbatchfile}\n")
                        print(hmm_cmd)
                        hmm_cmds.append(hmm_cmd)
        if "param_variation" in str(data_dir):
            for sel_type in ["add", "dom", "rec"]:
                for hs in [100, 250, 500, 1000, 2000]:
                    if hs == 1000:
                        mem = large_mem_required
                    elif hs == 2000:
                        mem = huge_mem_required
                    else:
                        mem = small_mem_required
                    in_name = Path(f"data/{sel_type}_s025_g251_d25_data.csv")
                    out_name = str(Path(EM_dir, f"{sel_type}_s025_g251_d25_hs{hs}_EM"))
                    sbatchfile = Path(qsub_dir, f"{sel_type}_s025_g251_d25_hs{hs}_EM.sbatch")
                    sbatchOutFile = Path(qsub_dir, f"{sel_type}_s025_251_d25_hs{hs}_EM.sout")
                    sbatchErrFile = Path(qsub_dir, f"{sel_type}_s025_g251_d25_hs{hs}_EM.serr")
                    hmm_cmd = f"emsel {in_name} {out_name} --time_after_zero --full_output -maf 0 --min_sample_density 0 --num_cores {num_cores} --selection_modes {sel_type} --no_neutral -hs {hs} --progressbar"
                    if (Path(out_name).with_suffix(".csv").is_file() and Path(out_name).with_suffix(
                            ".csv").stat().st_size > 0):
                        print(f"File already exists: {out_name}! continuing.")
                        continue
                    if hmm_cmd not in hmm_cmds:
                        # write a script for the cluster
                        with open(sbatchfile, 'w') as qsubScriptFile:
                            qsubScriptFile.write("#!/bin/bash\n" + \
                                                 f"#SBATCH -J {out_name[:10]}\n" + \
                                                 f"#SBATCH --time={max_run_hours}:00:00\n" + \
                                                 f"#SBATCH --cpus-per-task={num_cores}\n" + \
                                                 f"#SBATCH --mem={mem}gb\n" + \
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