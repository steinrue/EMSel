import pathlib

###### MODIFY
genodata_type = "capture_only"
max_run_hours = 10
num_cores = 30
hmm_script = pathlib.Path('../../run_emsel.py')
output_dir = pathlib.Path('output')
data_dir = pathlib.Path('data')

#this directory must be made beforehand
qsub_dir = pathlib.Path('qsubs')

###### DO NOT MODIFY

chroms = range(1, 23)
new_ugs = [["neutral", "add"], ["dom"], ["het"], ["rec"], ["full"]]
total_ugs = len(new_ugs)

def writeQsubs():
    hmm_cmds = []

    # get a meta file going
    script_file = pathlib.Path("meta_real_EM.sh")
    with open(script_file, "w") as file:
        for chrom in chroms:
            for ug_i, update_group in enumerate(new_ugs):
                in_name = pathlib.Path(data_dir, f"GB_v54.1_{genodata_type}_c{chrom}.vcf")
                out_name = pathlib.Path(output_dir, f"GB_v54.1_{genodata_type}_c{chrom}_EM_{ug_i+1}o{total_ugs}.pkl")
                sbatchfile = pathlib.Path(qsub_dir, f"GB_v54.1_{genodata_type}_c{chrom}_IM_{ug_i+1}o{total_ugs}.sbatch")
                sbatchOutFile = pathlib.Path(qsub_dir, f"GB_v54.1_{genodata_type}_c{chrom}_EM_{ug_i+1}o{total_ugs}.sout")
                sbatchErrFile = pathlib.Path(qsub_dir, f"GB_v54.1_{genodata_type}_c{chrom}_EM_{ug_i+1}o{total_ugs}.serr")
                # put together a nice command?
                hmm_cmd = f"python run_EMSel.py {in_name} {out_name} time_before_present --info_file {data_dir}/GB_v54.1_capture_only_inds.table --info_cols Genetic_ID Date_mean -ytg 28.1 --save_csv --full_output --num_cores {num_cores} --update_types {' '.join(u_i for u_i in update_group)} --no_neutral"

                if out_name.is_file() and out_name.stat().st_size > 0:
                    print(f"already exists: {out_name.name} {out_name.stat().st_size}")
                    continue
                if hmm_cmd not in hmm_cmds:
                    # write a script for the cluster
                    with open(sbatchfile, 'w') as qsubScriptFile:
                        qsubScriptFile.write("#!/bin/bash\n" + \
                                             f"#SBATCH -J {chrom}\n" + \
                                             f"#SBATCH --time={max_run_hours}:00:00\n" + \
                                             f"#SBATCH --cpus-per-task={num_cores}\n" + \
                                             "#SBATCH --mem=8gb\n" + \
                                             f"#SBATCH -o {sbatchOutFile}\n" + \
                                             f"#SBATCH -e {sbatchErrFile}\n" + \
                                             f"{hmm_cmd}\n")
                    # and add to meta file
                    file.write(f"sbatch {sbatchfile}\n")
                    hmm_cmds.append(hmm_cmd)


def main():
    writeQsubs()


if __name__ == "__main__":
    main()