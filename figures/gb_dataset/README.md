# GB aDNA Dataset

The scripts in this directory recreate Figures 12-14 of the main text and S1.R-AB, S1.AC(b), and S1.AD-AM of the Supplementary Material, as well as Table 1 of the main text and Tables S1.A-D of the Supplementary Material. This README begins with a set of commands that must be run prior to all figures, then  is organized by which (if any) additional non-plotting scripts must be run to generate a given set of figures. All commands to analyze and produce figures are described therein.

**Prior to running the scripts in this directory, change the current working directory to this folder.** Also, if you have not already, install the additional plotting packages required by running the command
```
pip install "emsel[plots] @ git+https://github.com/steinrue/EMSel"
```

## All figures

Each figure and table from this section of the manuscript and Supplementary Material requires a full analysis of the GB dataset described in Section 4.1 of the main text. UPDATE (11/3/25): we have added the full aggregated data file to Box, so that further analyses can be run without having to redo the full genome-wide scan! This file, which can be found [here](https://uchicago.box.com/s/mphspqmdwbg3tiesv0yuski1gp606ef7) (150 MB), contains the results of running the 8 steps below, as well as the results for the unconstrained EM (generated in the "Figures 14, S1.W, S1.AA, Table S1.D" section below). This file can be loaded into Python via the following snippet:
```
import pickle
complete_agg_data_path = **insert path to GB_v54.1_capture_only_agg_data.pkl file here**
with open(complete_agg_data_path, "rb") as file:
    data = pickle.load(file)
```

The `data` object will be a dictionary with the following keys:
- 'all_p' (dict): arrays of p-values computed for each SNP for each mode of selection in ["add", "dom", "het", "rec", "full"].
- 'all_s' (dict): arrays of selection coefficients estimated for each SNP for each mode of selection in ["add", "dom", "het", "rec", "full"]. Note that for the `full_s` key, corresponding to unconstrained EM, the array is two-dimensional, as both s\_1 and s\_2 are estimated.
- 'all_ll' (dict): arrays of log-likelihoods computed for each SNP for each mode of selection in ["add", "dom", "het", "rec", "full", "neutral"].
- 'all_loc_all_chrom' (array of ints, shape =(743417,)): genome-wide position for each SNP.
- 'all_loc_per_chrom' (array of ints, shape =(743417,)): position on a given chromosome for each SNP.
- 'all_chrom' (array of ints, shape =(743417,)): chromosome for each SNP.
- 'all_rsid' (array of strs, shape =(743417,)): rsid for each SNP.
- 'all_loc' (dict): contains, for each chromosome, the offset in indices and positions of SNPs in our dataset from 0 (i.e. the first SNP/beginning of chromosome 1), as well as arrays of the positions of all SNPs on a given chromosome for each SNP. the `all_loc_per_chrom` key in the main dictionary is an aggregation of each of these sub-arrays.
- 'all_missingness' (array of floats, shape =(743417,)): fraction of samples missing for each SNP.
- 'all_means' (array of floats, shape =(743417,)): initial frequency estimated under neutrality for each SNP.
- 'all_ref_allele' (array of strs, shape =(743417,)): reference allele at each SNP.
- 'all_alt_allele' (array of strs, shape =(743417,)): alternate allele at each SNP.
- 'maf_thresh' (float): minor allele frequency threshold used to filter the SNPs.
- 'missing_thresh' (float): missingness threshold used to filter the SNPs.

If you would like to recreate this file from scratch, proceed via the following:

1. Create the subfolders `data`, `EM`, `EM/Ne`, and `output`, and `qsubs` within this subdirectory.
2. Follow all instructions in the [extract_vcfs/](extract_vcfs/) subfolder. You should now have 44 .vcf files labelled capture_only_c{chr} and capture_SG_c{chr} for chr = (1,2,...,22) (the cX and cY files can be ignored from this point forward), as well as several additional .table files, in the extract_vcfs/extracted subfolder. In addition, step 4 of this process will have created Figures S1.R-S.
3. Move the contents of the `extract_vcfs/extracted` subfolder to the `data` subfolder created in step 1.
4. Run `python SLURM_Ne_final.py` with `EM_dir = Path('EM/Ne')`, `data_dir = Path('data')`, and `genodata_type = "capture_only"` at the beginning of the script. Modify other parameters to your liking. This should produce 462 jobs to submit to the cluster, corresponding to running a grid of 21 Ne values on each of the 22 autosomes. Then run `sh meta_gb_EM.sh` to submit the jobs to the cluster.
5. Run `python combine_split_runs.py` with `EM_dir = 'EM/Ne'` at the beginning of the script, then run `python interpolate_Ne_grid.py` to get a value of Ne estimated from the dataset. You should get a value of 12202. Now that the value of Ne has been estimated, we can run the full EM-HMM procedure:
6. Run `python SLURM_example.py` with `EM_dir = Path('EM')`, `data_dir = Path('data')`, and `genodata_type = "capture_only"` at the beginning of the script. Modify the other parameters to your liking. This should produce 110 jobs to submit to the cluster. Then, run `sh meta_gb_EM.sh` to submit the jobs to the cluster. Alternatively, for each created VCF whose filename contains `capture_only`, run EMSel via the command `emsel data/{file_name}.vcf EM/GB_v54.1_capture_only_c{chr}_EM --time_before_present --info_file data/GB_v54.1_capture_only_inds.table --info_cols Genetic_ID Date_mean -ytg 28.1 --save_csv --full_output`. Step 7 assumes that you have 22 files in `EM` that are named `GB_v54.1_capture_only_c{chr}_EM.pkl` for chr = (1,2,...,22).
7. Run `python combine_split_runs.py`, with `EM_dir = "EM"` at the beginning of the script. 
8. Set the following parameters at the beginning of the script `aggregate_data.py` and run it using `python aggregate_data.py`:
```
data_dir = "data"
EM_dir = "EM"
output_dir = "output"
genodata_type = "capture_only"
classification_types = ["add", "dom", "het", "rec", "full"]
```

## Figures 12, 13A, S1.X-Z, S1.AD-AK(a), Table 1*, Tables S1.A-C*

Set the following parameters at the beginning of the script `gb_figures.py` and run it using `python gb_figures.py`:
```
data_dir = "data"
EM_dir = "EM"
output_dir = "output"
genodata_type = "capture_only"
classification_types = ["add", "dom", "rec", "het"]
```

Note that the asterisk after the tables denotes that the generated tables are lacking the correct data in the "Gene(s)" column and confidence intervals for the selection coefficient estimates. These two datapoints are manually added post-hoc, though the script to simulate and analyze the data for the bootstrapped intervals, `sim_bootstraps.py`, is provided.

## Figures 13B, S1.AL-AM(a)

Set the following parameters at the beginning of the script `mathieson_manhattans.py` and run it using `python mathieson_manhattans.py`:
```
output_dir = "output"
genodata_type = "capture_only"
```

## Figures 13C+D, S1.AD-AM(b)

Set the following parameters at the beginning of the script `plot_binned_trajectories.py` and run it using `python plot_binned_trajectories.py`:
```
data_dir = "data"
EM_dir = "EM"
output_dir = "output"
genodata_type = "capture_only"
classification_types = ["add", "dom", "rec", "het"]
```

## Figures 14, S1.W, S1.AA, Table S1.D

To generate these figures, the p-values for the multiple alternatives procedure must be added to the aggregated data .pkl. To do this, set the following parameters at the beginning of the script `add_full_agg.py` and run it using `python add_full_agg.py`:
```
data_dir = "data"
EM_dir = "EM"
output_dir = "output"
genodata_type = "capture_only"
```

This will add the necessary sub-dictionaries to the `agg_data.pkl` file. Then, re-run `python gb_figures.py`, changing the classification_types line to `classification_types = ["full"]` at the beginning of the script.

## Figures S1.T-U
Set the following parameters at the beginning of the script `plot_means_and_missingness.py` and run it using `python plot_means_and_missingness.py`:
```
output_dir = "output"
genodata_type = "capture_only"
```

## Figure S1.V

Repeat steps 4-5 of the "All figures" pipeline, replacing `capture_only` with `capture_SG` everywhere it appears, and adding the option `--selection_modes neutral add` to the `emsel` command if you are not using the SLURM script. Then, rerun `python combine_split_runs.py`, `python aggregate_data.py` and `python gb_figures.py`, the latter two with `genodata_type = capture_SG` and `classification_types = ["add"]` substituted for their respective lines in the parameters at the top of the script.

## Figure S1.AB

Set the following parameters at the beginning of the script `plot_genomewide_correlations.py` and run it using `python plot_genomewide_correlations.py`:
```
output_dir = "output"
genodata_type = "capture_only"
```

## Figure S1.AC(b)

Run `python permute_gb_data.py` with the following parameters at the top of the script:
```
data_dir = "data"
output_dir = "output"
genodata_type = "capture_only"
MAF_filter = .05
min_sample_filter = .1
```
Then, for each created file (there will be 23 of them, all containing "permuted"), run `python SLURM_example.py` with the same parameters as the "All figures" section, submit the jobs using `sh meta_gb_em.sh`, and combine the runs afterwards with `python combine_split_runs.py` (with `EM_dir = "EM"`). Or, run EMSel via the command `emsel data/GB_v54.1_capture_only_c{chr}_permuted.csv EM/GB_v54.1_capture_only_c{chr}_permuted_EM --time_after_zero --full_output --selection_modes neutral add` for all chromosomes.

Lastly (or if you have already generated the permutation .pkl files), set the following parameters at the beginning of the script `plot_gb_permutations.py` and run it using `python plot_gb_permutations.py`:
```
EM_dir = "EM"
output_dir = "output"
genodata_type = "capture_only"
```

## Bootstrapping confidence intervals

To obtain bootstrapped confidence intervals for the estimated selection coefficients in Table 1 and Tables S1.A-C:
1. Create the subfolders `data/bootstrap`, `EM/bootstrap`, and `output/bootstrap`.
2. Set the following parameters at the beginning of the script `sim_gb_bootstraps.py` and run it using `python sim_gb_bootstraps.py`:
```
data_dir = "data/bootstrap"
output_dir = "output/boostrap"
genodata_type = "capture_only"
```
3. Run `python SLURM_example.py` with `EM_dir = Path('EM/bootstrap')` and `data_dir = Path('data/bootstrap')` at the beginning of the script. Submit the jobs by running `sh meta_gb_EM.sh`. Alternatively, For each .csv file created, run `emsel data/bootstrap/{file_name}.csv EM/{file_name}_EM --time_after_zero --full_output`
4. Run `python combine_split_runs.py` with `EM_dir = "EM/bootstrap"` at the beginning of the script.
5. Set the following parameters at the beginning of the script `compute_bootstraps.py` and run it using `python compute_bootstraps.py`:
  ```
data_dir = "data/bootstrap"
EM_dir = "EM/bootstrap"
output_dir = "output/bootstrap"
genodata_type = "capture_only"
``` 
This generates a set of boxplots, one for each selection mode, where the bias-corrected mean and confidence intervals for each significant SNP can be read off of its respective boxplot.
