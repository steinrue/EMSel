from pathlib import Path
import numpy as np

###### MODIFY

subdirectory = "ibdne"
data_dir = Path(f'data/{subdirectory}/boxplots')

qsub_dir = Path('qsubs')
max_run_hours = 1
num_cores = 1
huge_mem_required = 50
large_mem_required = 20
small_mem_required = 10

Ne = 9715

###### DO NOT MODIFY

if subdirectory == "ibdne":
    cmd_end = "--data_matched data/GB_v54.1_capture_only_means.txt data/GB_v54.1_capture_only_missingness.txt data/GB_v54.1_capture_only_sample_sizes.table --vary_Ne data/ibdne/ibdne_raw.txt"
    subdir_str = "ibdne"
elif subdirectory == "real_matched":
    cmd_end = f"--data_matched data/GB_v54.1_capture_only_means.txt data/GB_v54.1_capture_only_missingness.txt data/GB_v54.1_capture_only_sample_sizes.table -Ne {Ne}"
    subdir_str = "realmatch"

suffix = "unif"
cmd_suffix = " --ic_update_type fixed"

def writeQsubs():
    hmm_cmds = []
    # get a meta file going
    script_file = Path("meta_sim_EM.sh")
    with open(script_file, "w") as file:
                for i in range(25):
                    sbatchfile = Path(qsub_dir, f"sim_seed{100+i}.sbatch")
                    sbatchOutFile = Path(qsub_dir, f"sim_seed{100 + i}.sout")
                    sbatchErrFile = Path(qsub_dir, f"sim_seed{100 + i}.serr")
                    hmm_cmd = f"emsel-sim {data_dir} -s .005 --sel_types neutral --seed {100+i} -n 10000 --save_plots --suffix {subdir_str}_seed{100+i} {cmd_end}"
                    if hmm_cmd not in hmm_cmds:
                        # write a script for the cluster
                        with open(sbatchfile, 'w') as qsubScriptFile:
                            qsubScriptFile.write("#!/bin/bash\n" + \
                                                 f"#SBATCH -J Ne_{subdir_str}\n" + \
                                                 f"#SBATCH --time={max_run_hours}:00:00\n" + \
                                                 f"#SBATCH --cpus-per-task={num_cores}\n" + \
                                                 f"#SBATCH --mem={small_mem_required}gb\n" + \
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