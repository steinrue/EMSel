This directory contains scripts to extract the subset of the individuals in the GB from the AADR used for the analysis detailed in the manuscript. This subset comprises all individuals with:
- `political entity` given as `United Kingdom`.
- `Date` being between 4450 before present and 0 before present.
- `ASSESSMENT` given as `PASS`.
- No relatives in the dataset.
- Individuals from coastal islands removed manually.

We extract two datasets: One where `Data Source` is 1240K capture and shotgun, and one where `Data Source` is 1240K capture only.

The scripts use the the 1240K data from the AADR to extract genotypes (required) and the HO data from the AADR to produce PCA plots (optional). The scripts are configured to work with version v54.1 of the AADR.

The scripts are configured to expect the AADR to be located in a subdirectory named `AADR/`. To this end, download both the HO and 1240K data from the AADR v54.1 at https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/index_v54.1.html . Then place the respective files in the subdirectory as follows:
```
AADR/v54.1_HO/v54.1_HO_public.anno
AADR/v54.1_HO/v54.1_HO_public.geno
AADR/v54.1_HO/v54.1_HO_public.ind
AADR/v54.1_HO/v54.1_HO_public.snp
AADR/v54.1_1240K/v54.1_1240K_public.snp
AADR/v54.1_1240K/v54.1_1240K_public.anno
AADR/v54.1_1240K/v54.1_1240K_public.geno
AADR/v54.1_1240K/v54.1_1240K_public.ind
```
If the AADR is already downloaded on your device and stored at a different location, you can alternatively set up a symbolic link using: `ln -s <path/to/AADR> AADR`

Lastly, the scripts create the output files in a subdirectory named `extracted/`, which must exists and be empty in order for the scripts to run successfully. To create the empty subdirectory, execute the command: `mkdir extracted`

The scripts require that the programs `python` and `vcftools` are available to be called from the command line. Optionally, if you want to perform PCA, the tools `convertf` and `smartpca` from the `eigensoft` software suite need to be available. In addition, the python packages `pandas`, `numpy`, `collections`, `pathlib2`, `subprocess`, `matplotlib` (only for plotting), and `cartopy` (only for map-plot) need to be available.

All scripts define several global variables at the beginning of the files, that can be used to customize input and output.

#### 1. Extract annotation (required)

Exectuing:
```
python step1_subset_gb_anno.py
```
extracts the annotation for the subset of individuals to be analyzed.

#### 2. Prepare PCAs (optional)

*Steps 2, 3, and 4 are optional and only required if a PCA should be performed.*

The command:
```
python step2_eigensoft_convertf_subset.py
```
prepares all necessary input files for the PCA

#### 3. Run PCAs (optional)

*Steps 2, 3, and 4 are optional and only required if a PCA should be performed.*

The command:
```
python step3_eigensoft_smartpca_subset.py
```
performs the PCA.

Note that the script is configured to run all different PCAs sequentially, which can take while. If a compute cluster with the `slurm` queueing system is available, this can be parallelized. To this end, set the global variable `RUN_COMMANDS` in the script `step3_eigensoft_smartpca_subset.py` to `False`. Then the script will not execute the commands, but only print them to the console. Store these command in a file, for example using `python step3_eigensoft_smartpca_subset.py >SMARTPCA_COMMANDS.txt`. Then executing `python slurmify.py SMARTPCA_COMMANDS.txt` prepares scripts to be submitted, and by executing `sh metasubmit.sh` they will be send to the compute cluster.

#### 4. Plot PCAs (optional)

*Steps 2, 3, and 4 are optional and only required if a PCA should be performed.*

Running:
```
python step4_plot_pcas.py
```
will produce several PCA plots:
```
extracted/GB_v54.1_capture_only_broad_pca_noshrinkage.pdf
extracted/GB_v54.1_capture_only_broad_pca_shrinkage.pdf
extracted/GB_v54.1_capture_only_europe_pca_noshrinkage.pdf
extracted/GB_v54.1_capture_only_europe_pca_shrinkage.pdf
extracted/GB_v54.1_capture_only_gbr_ceu_pca_noshrinkage.pdf
extracted/GB_v54.1_capture_only_gbr_ceu_pca_shrinkage.pdf
extracted/GB_v54.1_capture_SG_broad_pca_noshrinkage.pdf
extracted/GB_v54.1_capture_SG_broad_pca_shrinkage.pdf
extracted/GB_v54.1_capture_SG_europe_pca_noshrinkage.pdf
extracted/GB_v54.1_capture_SG_europe_pca_shrinkage.pdf
extracted/GB_v54.1_capture_SG_gbr_ceu_pca_noshrinkage.pdf
extracted/GB_v54.1_capture_SG_gbr_ceu_pca_shrinkage.pdf
```
which are the results of the different PCAs. The files with `capture_only` result from using only the individuals genotyped on the 1240K panel, whereas `capture_SG` are the individuals genotype on the 1240K panel, as well as the shotgun individuals. `broad`, `europe`, and `gbr_ceu` refer to the 1000G reference populations used, which are {`LWK`, `YRI`, `JPT`, `CHB`, `CEU`, `GBR`, `FIN`, `TSI`, `IBS`}, {`CEU`, `GBR`, `FIN`, `TSI`, `IBS`}, and {`CEU`, `GBR`}, repsectively. Lastly, the analysis was performed using `shrinkage` or `noshrinkage`.

#### 5. Filter PCA outliers (required)

The command:
```
python step5_filter_pca_outlier_manual.py
```
will remove outliers identified by inspecting the PCAs. Currently, this actually does not remove any individuals, but it renames necessary files for downstream analyses.

#### 6. Plot maps indiciating sampling locations (optional)

The command:
```
python step6_plot_map.py
```
produces the plots
```
extracted/GB_v54.1_capture_only_map.pdf
extracted/GB_v54.1_capture_SG_map.pdf
```
that indicate the sampling locations of the individuals in the **capture only** or the **capture + shotgun** dataset.

#### 7. Plot sampling times and sizes (optional)

The command:
```
python step7_plot_sample_sizes.py
```
produces the plots
```
extracted/GB_v54.1_capture_only_sample_sizes.pdf
extracted/GB_v54.1_capture_SG_sample_sizes.pdf
```
that indicate the sampling times and sizes in the **capture only** or the **capture + shotgun** dataset.

#### 8. Extract genotypes for capture and shotgun individuals (required)

Executing:
```
python step8_extract_vcfs.py
```
extracts the genotypes of all individuals (**1240K capture** and **shotgun**) and stores them as VCFs. Specificially, this script produces a file
```
extracted/GB_v54.1_capture_SG_inds.table
```
that contains the annotations for all extracted individuals (**1240K capture** and **shotgun**), most importanly, their `Genetic_ID` to cross reference with the VCF and their sampling time `Date_mean`. In addition, the script generates the files
```
extracted/GB_v54.1_capture_SG_c*.vcf
```
with the genotype data, split up by chromosome. Note that the VCF files potentially have haploid genotypes for some and diploid genotypes for other individuals.

#### 9. Subset to capture only individuals (required)

Executing:
```
python step9_extract_capture_only.py
```
subsets the data from the previous step to only the individuals with **1240K capture** genotypes. Specificially, this script produces the file
```
extracted/GB_v54.1_capture_only_inds.table
```
that contains the annotations for the **1240K capture** individuals. In addition, the script generates the files
```
extracted/GB_v54.1_capture_only_c*.vcf
```
with the genotype data, split up by chromosome. Note that the VCF files potentially have haploid genotypes for some and diploid genotypes for other individuals.
