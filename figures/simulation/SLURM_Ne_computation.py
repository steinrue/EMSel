from pathlib import Path
import numpy as np

###### MODIFY
EM_dir = Path('EM/real_matched/')
data_dir = Path('data/real_matched/')

qsub_dir = Path('qsubs')
max_run_hours = 5
num_cores = 10
huge_mem_required = 50
large_mem_required = 20
small_mem_required = 10


###### DO NOT MODIFY

Nes = np.geomspace(500, 500000, num=21)
new_ugs = [["neutral", "add"], ["dom"], ["het"], ["rec"], ["full"]]
total_ugs = len(new_ugs)

suffix = "unif"
cmd_suffix = " --ic_update_type fixed"

def writeQsubs():
    hmm_cmds = []
    # get a meta file going
    script_file = Path("meta_sim_EM.sh")
    with open(script_file, "w") as file:
        for fpath in data_dir.iterdir():
            if fpath.suffix == ".csv":
                if "neutral" not in fpath.name:
                    continue
                for Ne in Nes:
                    out_name = str(Path(EM_dir, f"{fpath.name.rpartition('_')[0]}_HMMNe{int(Ne)}_cond_{suffix}_EM"))
                    sbatchfile = Path(qsub_dir, f"{fpath.name.rpartition('_')[0]}_HMMNe{int(Ne)}_cond_{suffix}_EM.sbatch")
                    sbatchOutFile = Path(qsub_dir, f"{fpath.name.rpartition('_')[0]}_HMMNe{int(Ne)}_cond_{suffix}_EM.sout")
                    sbatchErrFile = Path(qsub_dir, f"{fpath.name.rpartition('_')[0]}_HMMNe{int(Ne)}_cond_{suffix}_EM.serr")
                    hmm_cmd = f"emsel {fpath} {out_name} --time_after_zero --full_output --num_cores {num_cores} --selection_modes neutral --progressbar --compute_cond -Ne {int(Ne)}{cmd_suffix}"

                    if Path(out_name+'.csv').is_file() and Path(out_name+".csv").stat().st_size > 0:
                        print(f"File already exists: {out_name}! continuing.")
                        continue

                    if hmm_cmd not in hmm_cmds:
                        # write a script for the cluster
                        with open(sbatchfile, 'w') as qsubScriptFile:
                            qsubScriptFile.write("#!/bin/bash\n" + \
                                                 f"#SBATCH -J {fpath.name[:5]}\n" + \
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