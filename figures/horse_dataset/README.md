# ASIP locus analyses

The scripts in this directory recreate Figure 15 and Table 2 of the main text, as well as the accompanying results regarding mode inference.

**Prior to running the scripts in this directory, change the current working directory to this folder.** Also, if you have not already, install the additional plotting packages required by running the command
```
pip install "emsel[plots] @ git+https://github.com/steinrue/EMSel"
```

## Figure 15A

Run `python plot_horse_trajectory.py`. 

## Figures 15B+C, Table 2

Proceed via the following:

1. Create the subfolders `data`, `EM`, `output`, and `qsubs` within this subdirectory.
2. Move the files `horse_data_full.csv`, `horse_data_colon3.csv`, and `full_sample_sizes_8.table` from [sample_datasets/](../../sample_datasets/) to the `data` subfolder of this directory.
3. Run `python SLURM_example.py`, followed by `sh meta_sim_EM.sh`. Alternatively, for each .csv file now in the `data` directory, run `emsel {data_path} {output_path} --time_after_zero --full_output -Ne 16000 --selection_modes all --progressbar`.
4. Next, run `python extract_means.py` to obtain the initial allele frequency estimates for the ASIP locus under neutrality and the unconstrained EM.
5. Run `python horse_final_simulation.py` to generate the simulated replicates for the empirical distribution of the D_2 statistic.
6. Re-run `python SLURM_example.py` followed by `sh meta_sim_EM.sh` to analyze these simulations.
7. Lastly, run `python horse_final_analysis.py` to generate Figures 15B+C and Table 2.
