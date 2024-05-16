# GB aDNA Dataset

The scripts in this directory recreate Figures 12-14 of the main text and S.12-S21, S.22B, S.23-S.32 of the Supplementary Material, as well as Table 1 of the main text and Tables S.1-S.3 of the Supplementary Material. This README begins with a set of commands that must be run prior to all figures, then  is organized by which (if any) additional non-plotting scripts must be run to generate a given set of figures. All commands to analyze and produce figures are described therein.

**Prior to running the scripts in this directory, change the current working directory to this folder.** Also, if you have not already, install the additional plotting packages required by running the command
```
pip install "emsel[plots] @ git+https://github.com/steinrue/EMSel"
```

## All figures

Each figure and table from this section of the manuscript and Supplementary Material requires a full analysis of the GB dataset described in Section 4.1 of the main text. To do so, proceed via the following:

1. Create the subfolders `data`, `EM` and `output` within this subdirectory.
2. Follow all instructions in the [extract_vcfs/](extract_vcfs/) subfolder. You should now have 48 .vcf files, as well as several additional .table files, in the extract_vcfs/extracted subfolder. In addition, step 4 of this process will have created Figures S.12-S.13.
3. Move the contents of the `extract_vcfs/extracted` subfolder to the `data` subfolder created in step 1.
4. For each created VCF whose filename contains `capture_only`, run EMSel via the command `emsel data/{file_name}.vcf EM/GB_v54.1_capture_only_c{chr}_EM --time_before_present --info_file data/GB_v54.1_capture_only_inds.table --info_cols Genetic_ID Date_mean -ytg 28.1 --save_csv --full_output`. Due to the size of the VCFs, this step is quite computationally expensive - we recommend running this in parallel on a cluster (using the -nc flag) and splitting the runs by selection mode. The scripts `SLURM_example.py` and `combine_split_runs.py` provide a template for submitting scripts that are parallelized and split by selection mode to a cluster and combining the results into a single file afterwards, respectively. Step 5 assumes that you have 22 files in `EM` that are named `GB_v54.1_capture_only_c{chr}_EM.pkl` for chr = (1,2,...,22).
5. Set the following parameters at the beginning of the script `aggregate_data.py` and run it using `python aggregate_data.py`:
```
data_dir = "data"
EM_dir = "EM"
output_dir = "output"
genodata_type = "capture_only"
```
This will also generate the `GB_v54.1_capture_only_means.txt` and `GB_v54.1_capture_only_missingness.txt` files needed for the data-matched simulations (see "Figures 9-11" in the [figures/simulation](../simulation) README). The data-matched simulations, in turn, are needed to analyze the unconstrained EM and recreate Figure 14.

## Figures 12, 13A, S.17-20, (S.23-30)A, Table 1*, Tables S.1-S.3*

Set the following parameters at the beginning of the script `gb_figures.py` and run it using `python gb_figures.py`:
```
data_dir = "data"
EM_dir = "EM"
output_dir = "output"
genodata_type = "capture_only"
classification_types = ["add", "dom", "rec", "het"]
```

Note that the asterisk after the tables denotes that the generated tables are lacking the correct data in the "Gene(s)" column and confidence intervals for the selection coefficient estimates. These two datapoints are manually added post-hoc, though the script to simulate and analyze the data for the bootstrapped intervals, `sim_bootstraps.py`, is provided.

## Figures 13B, (S.31-32)B

Set the following parameters at the beginning of the script `mathieson_manhattans.py` and run it using `python mathieson_manhattans.py`:
```
output_dir = "output"
genodata_type = "capture_only"
```

## Figures 13C+D, (S.23-30)B

Set the following parameters at the beginning of the script `plot_binned_trajectories.py` and run it using `python plot_binned_trajectories.py`:
```
data_dir = "data"
EM_dir = "EM"
output_dir = "output"
genodata_type = "capture_only"
classification_types = ["add", "dom", "rec", "het"]
```

## Figure 14

To generate this figure, the unconstrained EM must be classified. For this, the `gengamma_params.pkl` file is needed. To generate this file, either:
1. Run the scripts in the "Figures 9-11" section of the [figures/simulation](../simulation) folder (everything before the "Figure 9A" header), followed by the `deltall_qqs_and_confusiontables.py` script as detailed in the "Figure 9D+10" section.
2. Move the provided `gengamma_params.pkl` file from the [sample_datasets](../../sample_datasets) folder into `data`.

Then, set the following parameters at the beginning of the script `add_full_agg.py` and run it using `python add_full_agg.py`:
```
data_dir = "data"
EM_dir = "EM"
output_dir = "output"
genodata_type = "capture_only"
```

This will add the necessary sub-dictionaries to the `agg_data.pkl` file. Then, re-run `python gb_figures.py`, changing the classification_types line to `classification_types = ["full"]` at the beginning of the script.

## Figures S.14+S.15

Set the following parameters at the beginning of the script `plot_means_and_missingness.py` and run it using `python plot_means_and_missingness.py`:
```
output_dir = "output"
genodata_type = "capture_only"
```

## Figure S.16

Repeat steps 4-5 of the "All figures" pipeline, replacing `capture_only` with `capture_SG` everywhere it appears, and adding the option `--selection_modes neutral add` to the `emsel` command. Then, rerun `python gb_figures.py` with `genodata_type = capture_SG` and `classification_types = ["add"]` substituted for their respective lines in the parameters at the top of the script.

## Figure S.21

Set the following parameters at the beginning of the script `plot_genomewide_correlations.py` and run it using `python plot_genomewide_correlations.py`:
```
output_dir = "output"
genodata_type = "capture_only"
```

## Figure S.22B

First, set the following parameters at the beginning of the script `permute_gb_data.py` and run it using `python permute_gb_data.py`:
```
data_dir = "data"
output_dir = "output"
genodata_type = "capture_only"
MAF_filter = .05
min_sample_filter = .1
```

Then, for each created file (there will be 22 of them, all containing "permuted"), run EMSel via the command `emsel data/GB_v54.1_capture_only_c{chr}_permuted.csv EM/GB_v54.1_capture_only_c{chr}_permuted_EM --time_after_zero --full_output --selection_modes neutral add`

Lastly, set the following parameters at the beginning of the script `plot_gb_permutations.py` and run it using `python plot_gb_permutations.py`:
```
EM_dir = "EM"
output_dir = "output"
genodata_type = "capture_only"
```

## Bootstrapping confidence intervals

To obtain bootstrapped confidence intervals for the estimated selection coefficients in Table 1 and Tables S.1-S.3:
1. Create the subfolders `data/bootstrap` and `EM/bootstrap`.
2. Set the following parameters at the beginning of the script `sim_gb_bootstraps.py` and run it using `python sim_gb_bootstraps.py`:
```
data_dir = "data/bootstrap"
output_dir = "output"
genodata_type = "capture_only"
```
3. For each .csv file created, run `emsel data/bootstrap/{file_name}.csv EM/{file_name}_EM --time_after_zero --full_output`
4. Set the following parameters at the beginning of the script `compute_bootstraps.py` and run it using `python compute_bootstraps.py`:
  ```
data_dir = "data/bootstrap"
EM_dir = "EM/bootstrap
output_dir = "output"
genodata_type = "capture_only"
``` 
This generates a set of boxplots, one for each selection mode, where the bias-corrected mean and confidence intervals for each significant SNP can be read off of its respective boxplot.
