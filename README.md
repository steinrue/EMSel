# EMSel
Code accompanying Fine & Steinrücken (2024): A novel expectation-maximization approach to infer general diploid selection from time-series genetic data. (https://doi.org/10.1101/2024.05.10.593575)

We provide the command-line tool `emsel` to analyze time-series allele frequency datasets under multiple modes of selection (additive, dominant, recessive, over/underdominance, general diploid). In addition, `emsel-sim` can be used to simulate such data.

The subdirectory [figures/](figures/) contains all code to compile the data and reproduce the figures from the paper.

## Installation

To install EMSel, run the following `pip` command on the unix command-line:
```
pip install "emsel @ git+https://github.com/steinrue/EMSel"
```
This will install all necessary packages to run the `emsel` and `emsel-sim` command-line tools. 

To additionally install packages needed to run the scripts in the [figures/](figures/) subdirectory, use the command
```
pip install "emsel[plots] @ git+https://github.com/steinrue/EMSel"
```
Also, note that all example scripts in this README assume that `sample_datasets` is a subfolder of the current working directory.


## Running EMSel

To run EMSel, you must have either a CSV file or a VCF file. 
### Using EMSel with CSVs

The CSV should be formatted as the following:
- each locus/replicate should be in its own row.
- each row should contain 3N non-zero values, consisting of N samples each formatted as the triple (sampling time, number of derived alleles, number of total samples).
- samples of size zero are allowed, and are interpreted identically to a gap between sampling times.
- rows do not need to have the same length.
- sampling times can be expressed in years or generations: see the use of the -ytg flag below.

An example of a properly-formatted CSV is available in `sample_datasets/add_s025_g251_d25_data.csv`

### Using EMSEL with VCFs

Any VCF that can be read by `scikit-allele`, it can be used with EMSel. Using EMSel with a VCF additionally requires a file containing the same strings as the `samples` key in the VCF and a corresponding sampling time (in years or generations, see the -ytg flag below) for each sample.

### Minimal example and output

A minimal sample call to EMSel with a CSV:
`emsel sample_datasets/add_s025_g251_d25_data.csv output_EM --time_after_zero --progressbar`


A minimal sample call with a VCF:
`emsel sample_datasets/GB_c22.vcf output_EM --info_file sample_datasets/GB_individuals.table --info_cols Genetic_ID Date_mean --time_before_present --progressbar`

Both of these will create the file `output_EM.csv` containing a simple table of the results of running EMSel in all available modes of selection. The table is formatted as a indexed and column-labelled csv with one row for each replicate with (3M+2) columns per row, where M equals the number of non-neutral selection modes analyzed under. The first column are the index of each row within the unfiltered dataset (see -maf and --min_sample_density for a description of the filters). The second column is the neutral log-likelihood. Each set of 3 subsequent columns is the tuple (log_likelihood, s_1, s_2) for each selection mode at termination of the algorithm.

A more complete output file will also be saved if the `--full_output` flag is used - see that section of the README for a description. 
 
### Command-line arguments

In addition to the required `input` and `output` paths, EMSel has one required argument:

```
--time_after_zero | --time_before_present
    Use of exactly one of these arguments is required to specify whether the sample dates provided in the CSV/VCF start
    at zero at the earliest time and count up towards the present (--time_after_zero) or start at zero at the present and
    count up moving backward in time (--time_before_present).
```

and the following optional arguments:

```
-ytg, --years_to_gen <float, default=1>
    Number of years per generation, used to convert a VCF or CSV to generations. If the sampling times in the sample
    file or the CSV are in generations, use the defualt of 1. Note that for the --save_csv flag, the CSV output
    will be in generations.


-maf, --min_allele_freq <float, default=0.05>
    Minor allele frequency (MAF) threshold to filter replicates by.
    Replicates with min(mean(ref_allele), mean(alt_allele)) < MAF are masked.


--min_sample_density <float, default=0.1>
    Minimum sample density to filter replicates by.
    Replicates with sum(num_samples) < (max_samples * min_sample_density_thresh) are masked. max_samples is
    calculated correctly (inferring haploids/diploids) from a VCF, and simply computed as max(sum(num_samples, axis=1))
    for CSVs, where the sum is taken across timepoints for each replicate.


--s_init <2 floats, default=0 0>
    (s1, s2) for initialization of the HMM and the inital "guess" for each replicate.


-sid, --starting_init_dist <str, default='uniform'>
    Starting initial distribution (pi_0 in HMM language) for initialization of the HMM. Options are: 
    - "uniform" - uniform, equiprobable prior (recommended/default)
    - "delta" - use the "--ic_dict p x", where 0 <= x <= 1 to set the initial condition to a delta function at p = x.
    - "beta" - use "--ic_dict beta_coef alpha", where alpha is a real number, to set the initial condition to
         a symmetric beta function beta(alpha, alpha).
    - "spikeandslab" - use "--ic_dict spike_frac x spike_loc y", where 0 <= x, y <= 1 to set the initial condition to a
         mixture of uniform with weight 1-x and a delta function at p = y with weight x.
    Unless --ic_update_type fixed is provided as an additional command line argument, this starting initial distribution
    will be re-estimated as part of the EM procedure.


--sid_dict <1+ arguments of the form 'str int'>
    Additional arguments if -sid is not "uniform" or "theta".


-Ne <int, default=10000>
    Effective population size.


--hidden_interp <str, default='chebyshev'>
    Whether to use Chebyshev nodes for spacing of hidden states (highly recommended) or linear
    (via --hidden_interp linear). Chebyshev nodes do not impact runtime and appear to significantly improve accuracy
    of selection coefficient estimation at high selection coefficients, especially for certain modes of selection.


--ic_update_type <str, default='beta'>
    Method of estimating the initial condition. Options are:
    - "beta" - estimate the parameters of a beta distribution. Output dictionary values involving
        the initial distribution will have shape (N,2)
    - "delta" - estimate the parameter of a delta distribution. Output dictionary values will have shape (N, 1)
    - "baumwelch" - use the standard Baum-Welch EM update rules to estimate the weights for all hidden states.
        Output dictionary values will have shape (N, Ns), where Ns = the number of hidden states.
    - "fixed" - do not update the estimated initial condition. Output dictionary values will have shape (N, 1)
    Any other string will raise a ValueError.


--selection_modes <1+ str, default='all'>
    Which modes of selection to analyze under. You can list as many of ["neutral", "add", "dom", "rec", "het", "full"]
    as you would like. Neutral is automatically run. The default, "all", is a shorthand for running all modes of
    selection.


--no_neutral
    Explicitly overrides the requirement that neutral be run for all simulation. Only to be used when parallelizing
    a simulation across mode of selection.


-nc, --num_cores <int, default=1>
    Number of cores to parallelize over. Joblib is used for parallelization, making EMSel easily parallizable on a
    computing cluster.


-hs, --hidden_states <int, default=500>
    Number of hidden states to use in the HMM. Computing time scales as O(hidden states^3), roughly.
    Accuracy is reduced below approximately 200-250 hidden states.


-t, --tol <float, default=1e-3>
    Stopping criterion - EMSel converges for a replicate if
    (log-likelihood at iteration k+1 - log-likelihood at iteration k) < tol.


-m, --maxiter <int, default=2000>
    Maximum number of iterations of EMSel before terminating.


--min_itercount <int, default=5>
    Minimum number of iterations for EMSel to run, even if the stopping criterion is met earlier. Removes irregularities
    in the distribution of log-likelihoods/p-values near 1. However, since p-values near 1 are typically unimportant,
    setting this to 0 to slightly speed up computation is reasonable.


--info_file <str>
    When input is a VCF, path to a sample file readable by pandas.read_csv containing a column of IDs matching the IDs
    in the VCF as well as a column of sampling times
    (in years or generations, use -ytg to normalize to generations if years).


--info_cols <2 strs>
    Column names, in (IDs, times) order, to extract from the sample file.


--save_csv
    If input is a VCF, saves a CSV of the same name containing the intermediate conversion of the VCF into
    (sampling time, total samples, derived alleles) triplets to speed up future runs. Note that the saved CSV will
    include conversion from years to generations.


--full_output
    Saves a full output file to the same location as the output_EM.csv file (with a .pkl suffix). The full output
    file is a nested dictionary-of-dictionaries containing the following keys, letting `N_init` and `N` be
    the number of input and filtered replicates, respectively, in the dataset:
    - `max_samples` (1,) - the maximum number of samples for a single replicate in the dataset.
    - `sample_mask` (N_init,) - a boolean mask, where True indicates that a given replicate is in the final dataset
         (i.e. has not been filtered out)
    - `sample_idxs (N,) - the indices of the filtered replicates with respect to the rows of the initial dataset.
    - `neutral_ll` (N,) - the log-likelihood for each replicate calculated under neutrality (s1 = s2 = 0).
    - `neutral_ic` (N, varies) - the estimated initial distribution parameters for each replicate calculated under
        neutrality. The second dimension depends on which initial distribution is used for calculation.
    - `neutral_itercount` (N,) - the number of iterations for convergence for each replicate under neutrality.
    - for each mode of selection analyzed under, a sub-dictionary with key `{update_type}_run`
        (e.g. `add_run` for additive selection), containing the following keys:
          - `s_final` (N, 2) - the maximum-likelihood estimate of the selection coefficients for each replicate.
          - `ll_final` (N,) - the maximum log-likelihood estimate for each replicate.
          - `ic_dist` (N, varies) - the estimated initial distribution parameters for each replicate. The second
                dimension depends on which initial distribution is used for calculation.
          - `itercount_hist` (N,) - the number of iterations for convergence for each replicate.
          - `exit_codes` (N,) - exit codes indicating the termination statuts of each replicate.
               See section "Exit Codes".

    Additionally, if the input file is a .vcf, several additional keys are added from the VCF file to facilitate
    additional data analysis (see the figures/gb_dataset folder), each of shape (N,):
    - `pos` - chromosomal position of each replicate. "variants/POS" from the VCF file.
    - `snp_ids` - "variants/ID" from the VCF file.
    - `ref_allele` and `alt_allele` - "variants/REF" and ["variants/ALT"][:, 0] from the VCf file.


--force [haploid|diploid]
    If the inputted VCF file contains only homozygous loci, use this flag to determine whether the genotypes are
    read as haploid (1 sample per locus per individual) or homozygous diploid (2 samples per locus per individual).
    If the inputted VCF file is not all homozygotes, this argument has no effect.


--progressbar
    Track the progress of EMSel with a tqdm progressbar.
```


### Exit Codes

The 'exit_codes' array contained in the full_output file has the following possible values for each replicate:
- 0: Replicated converged successfully.
- 1: Initial distribution estimation failed, causing immediate termination of the EM algorithm.
- 2: Invalid EM step proposed (s values too large). Initial distribution was set to fixed, so no improvements possible. Immediate termination of EM algorithm.
- 3: Decreasing likelihood for two or more iterations, outside of the first `min_itercount` iterations. Immediate termination of EM algorithm.
- 4: Replicate converged successfully, but at some point an invalid EM step was proposed. Initial distribution continued to be updated subsequently.
- 7: Continuing for `min_itercount` iterations caused the likelihood to decrease. Parameters at the iteration with maximum likelihood were returned.
- 8: Convergence successfully achieved before `min_itercount` iterations.
- 9: Maximum number of iterations exceeded without convergence.
- 12: Convergence successfully achieved before `min_itercount` iterations, with an invalid EM step at some iteration (combination of exit codes 4 and 8).

## Simulating data

`simulate_data.py` (used as `emsel-sim` from the command line) uses the discrete Wright-Fisher model to simulate allele frequencies under a given selection scenario (selection strength, mode of selection, initial condition, number of generations simulated, sampling scheme). 

### Output and examples

`emsel-sim` outputs three different file types:

1. `args_{suffix}.pkl` - contains a dictionary of the parameters used to run all sets of simulations. 

For each simulation condition, the following two files are outputted:

2. `{exp_name}_{suffix}data.csv` - contains the sampled allele frequencies and sampling times (see Using EMSel with CSVs for full formatting description).
3. `{exp_name}_{suffix}pd.bz2` - a bz2-zipped pickle file containing the simulation parameters for this particular simulation condition.

Sample calls to `simulate_data.py` for a non-data-matched set of simulations and a data-matched simulation are as follows:

`emsel-sim . -s .01 .1 -g 101 251 -ic .05 recip --suffix big_s`

and

`emsel-sim . -s .005 .05 .2 --sel_types neutral add rec --data_matched sample_datasets/GB_means.txt sample_datasets/GB_missingness.txt sample_datasets/GB_sample_sizes.table`

### Command-line arguments

There is one required argument - a path to an output directory. In addition, there are the following positional arguments:

```
-n, --num_sims <int, default=100>
Number of replicates simulated from each set of simulation parameters.


-s, --sel_coeffs <float or sequence of floats, default=.005, .01, .025, .05>
    Selection coefficients to simulate. For additive selection, this value is s_2 (the fitness of the derived homozygote),
    with s_1 = s_2 / 2. All other modes of selection are governed by a single parameter. Simulating under neutrality
    should be done by adding "neutral" to the --sel_types flag, not by including 0 in this flag.


--sel_types <str or sequence of strs, default=neutral add dom rec over under>
    Modes of selection to simulate under. Must be a (non-strict) subset of
    {"neutral", "add", "dom", "rec", "over", "under"}.


-g, --num_gens <int or sequence of ints, default=101 251 1001>
    Number of generations to simulate for. 


-ic, --init_conds <float, str, or sequence of floats and strs, default=.005, .25, "recip">
    Initial conditions to simulate under. Floats are interpreted as initializing each simulation from the same fixed
    allele frequency. Currently, the only non-float option is "recip", which samples from a distribution where the
    probability of initial frequency p is proportional to 1/p.



-ns, --num_samples <int, default=50>
    Number of haploid samples drawn at each sampling timepoint.



-st, --sampling_times <int, default=11>
    Number of equally-spaced timepoints to sample at. Samples are taken at the first generation,
    the last generation and (sampling_times - 2) points equally spaced between.



-Ne <int, default=10000>
    Effective population size for the simulations.


--data_matched <3 strs>
    Provide the path to a sample_means.txt file, a sample_missingness.txt file and a sample_sizes.table file to
    override the -ic, -g, -ns, and -st flags and simulate under the same initial distribution, number of generations,
    and sampling scheme as the inputted real data file. The sample_means.txt and sample_missingness.txt files can be
    obtained by running the aggregate_data.py script (in figures/gb_dataset) after running EMSel. sample_means.txt
    should be formatted as a MAF filter on the first line followed by one float per line representing estimated
    initial frequencies to draw from. sample_missingness.txt should be formatted as a missingness filter on the first
    line followed by one float per line representing sample missingness to draw from. The sample_sizes.table file
    can be obtained by running the pipeline in the figures/gb_dataset/extract_vcfs folder.
    Example files are also provided.


--no_small_s
    By default, the formula for updating allele frequencies in the Wright-Fisher model uses the small s approximation
    (p' = p + p(1-p)*((1-2p)s1 + p*s2)). Use this flag to indicate that the full formula,
    p' = p(1+s1(1-p)+s2*p)/(1+2*s1*p+s2*p^2-2s1*p**2), should be used instead.


--seed <int>
    Seed to simulate from. Simulations from identical seed values are identical.


--save_plots
    If used, one plot per simulation condition plotting the true allele frequency for each replicate as well as
    the mean allele frequency will be produced.


--suffix <str>
    Adds a suffix to each file to help distinguish different simulation runs.
```

When specifying multiple values for sel_coeffs, sel_types, num_gens, and init_conds, simulation is done for each combination on values (i.e. on itertools.product(sel_coeffs, sel_types, num_gens, init_conds)), for a total of len(sel_coeffs)\*len(sel_types)*(len(num_gens)*len(init_conds) sets of files. 
