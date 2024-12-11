from pathlib import Path

###### MODIFY
EM_dir = Path('EM')
data_dir = Path('data')

qsub_dir = Path('qsubs')
max_run_hours = 24
num_cores_default = 2
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
            if fpath.suffix == ".csv" and fpath.name[0] != ".":
                out_name = str(Path(EM_dir, f"{fpath.stem}_EM"))
                sbatchfile = Path(qsub_dir, f"{fpath.stem}_EM.sbatch")
                sbatchOutFile = Path(qsub_dir, f"{fpath.stem}_EM.sout")
                sbatchErrFile = Path(qsub_dir, f"{fpath.stem}_EM.serr")
                if "2tp" not in fpath.name:
                    continue
                if "perm" in fpath.name:
                    num_cores = 20
                    hmm_cmd = f"emsel {fpath} {out_name} --time_after_zero --full_output --num_cores {num_cores} -Ne 2500 --selection_modes neutral add full --progressbar"
                elif "Ne" in fpath.name:
                    strname = str(fpath.name)
                    num_cores = 30
                    Ne_loc = strname.index("Ne")
                    next_underscore = strname[Ne_loc+2:].index("_")
                    Ne_val = strname[Ne_loc+2:Ne_loc+2+next_underscore]
                    hmm_cmd = f"emsel {fpath} {out_name} --time_after_zero --full_output --num_cores {num_cores} -Ne {int(Ne_val)} --selection_modes all --progressbar -maf 0 --min_sample_density 0"
                else:
                    num_cores = num_cores_default
                    hmm_cmd = f"emsel {fpath} {out_name} --time_after_zero --full_output --num_cores {num_cores} -Ne 16000 --selection_modes all --progressbar"
                if (Path(out_name).with_suffix(".csv").is_file() and Path(out_name).with_suffix(
                        ".csv").stat().st_size > 0):
                    print(f"File already exists: {out_name}! continuing.")
                    continue

                if hmm_cmd not in hmm_cmds:
                    # write a script for the cluster
                    with open(sbatchfile, 'w') as qsubScriptFile:
                        qsubScriptFile.write("#!/bin/bash\n" + \
                                             f"#SBATCH -J {fpath.name[:5]}\n" + \
                                             f"#SBATCH --time={max_run_hours}:00:00\n" + \
                                             f"#SBATCH --cpus-per-task={num_cores}\n" + \
                                             f"#SBATCH --mem={large_mem_required}gb\n" + \
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