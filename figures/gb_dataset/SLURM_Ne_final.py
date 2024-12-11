from pathlib import Path
import numpy as np

###### MODIFY
EM_dir = Path('EM/Ne')
data_dir = Path('data')
genodata_type = "capture_only"

qsub_dir = Path('qsubs')
max_run_hours = 2
num_cores = 30
mem_required = 15


#this directory must be made beforehand


###### DO NOT MODIFY

chrom = 3

Nes = np.geomspace(500, 500000, num=21)
new_ugs = [["neutral", "add"], ["dom"], ["het"], ["rec"], ["full"]]
total_ugs = len(new_ugs)

def writeQsubs():
    hmm_cmds = []

    # get a meta file going
    script_file = Path("meta_gb_EM.sh")
    with open(script_file, "w") as file:
        for fpath in data_dir.iterdir():
            if genodata_type not in fpath.name or "binned" in fpath.name or "permuted" in fpath.name or fpath.suffix != ".csv":
                continue
            for Ne in Nes:
                out_name = str(Path(EM_dir, f"{fpath.stem}_Ne{int(Ne)}_cond_unif_EM"))
                sbatchfile = Path(qsub_dir, f"{fpath.stem}_Ne{int(Ne)}_cond_unif_EM.sbatch")
                sbatchOutFile = Path(qsub_dir, f"{fpath.stem}_Ne{int(Ne)}_cond_unif_EM.sout")
                sbatchErrFile = Path(qsub_dir, f"{fpath.stem}_Ne{int(Ne)}_cond_unif_EM.serr")
                # put together a nice command?
                hmm_cmd = f"emsel {fpath} {out_name} --num_cores {num_cores} --time_after_zero --full_output --selection_modes neutral --progressbar --compute_cond --ic_update_type fixed -Ne {int(Ne)}"

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