from pathlib import Path

###### MODIFY
max_run_hours = 10
num_cores = 30
EM_dir = Path('EM')
data_dir = Path('data')

#this directory must be made beforehand
qsub_dir = Path('qsubs')

###### DO NOT MODIFY
extra_fnames = ["Ne5000", "Ne20000","ns100_linear", "linear", "fixed_ic"]
extra_cmds = ["-Ne 5000", "-Ne 20000", "-hs 100 --hidden_interp linear --ic_update_type fixed",
              "--hidden_interp linear --ic_update_type fixed", "--ic_update_type fixed"]

def writeQsubs():
    hmm_cmds = []
    # get a meta file going
    script_file = Path("meta_sim_EM.sh")
    with open(script_file, "w") as file:
        for fpath in data_dir.iterdir():
            if fpath.suffix == ".csv":
                if "g251_d25" in fpath.name:
                    for e_i, extra_cmd in enumerate(extra_cmds):
                        mem=10
                        sel_type = fpath.name.split("_")[0]
                        if sel_type == "neutral":
                            sel_modes = "all"
                        elif sel_type == "add":
                            sel_modes = "neutral add"
                        elif sel_type == "over" or sel_type == "under":
                            sel_modes = "het"
                        else:
                            sel_modes = sel_type
                        out_name = str(Path(EM_dir, f"{fpath.name.rpartition('_')[0]}_{extra_fnames[e_i]}_EM"))
                        sbatchfile = Path(qsub_dir, f"{fpath.name.rpartition('_')[0]}_{extra_fnames[e_i]}_EM.sbatch")
                        sbatchOutFile = Path(qsub_dir, f"{fpath.name.rpartition('_')[0]}_{extra_fnames[e_i]}_EM.sout")
                        sbatchErrFile = Path(qsub_dir, f"{fpath.name.rpartition('_')[0]}_{extra_fnames[e_i]}_EM.serr")
                        hmm_cmd = f"emsel {fpath} {out_name} --time_after_zero --full_output -maf 0 --min_sample_density 0 --num_cores {num_cores} --selection_modes {sel_modes} --no_neutral {extra_cmd} --progressbar"
                        if (Path(out_name).with_suffix(".csv").is_file() and Path(out_name).with_suffix(".csv").stat().st_size > 0):
                            print(f"File already exists: {out_name}! continuing.")
                            continue

                        if hmm_cmd not in hmm_cmds:
                            # write a script for the cluster
                            with open(sbatchfile, 'w') as qsubScriptFile:
                                qsubScriptFile.write("#!/bin/bash\n" + \
                                                     f"#SBATCH -J {fpath.name[:10]}\n" + \
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