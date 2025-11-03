import numpy as np
from emsel.core import HMM
import pickle
from pathlib import Path
from tqdm import tqdm
import argparse
import allel
from pandas import read_csv
from emsel.emsel_util import vcf_to_useful_format, save_csv_output
from joblib import Parallel, delayed
from scipy.stats import beta

def run_one_s(iter_hmm, obs_counts, nts, sample_locs, loc, tol, max_iter, init_mean=False, min_init_val=1e-8, min_ic = 5, save_history=False):
    iter_hmm.update_internals_from_datum(obs_counts, nts, sample_locs)
    iter_hmm.s1 = iter_hmm.s1_init
    iter_hmm.s2 = iter_hmm.s2_init
    iter_hmm.a = iter_hmm.a_init
    if init_mean:
        mean_freq = np.sum(obs_counts)/np.sum(nts)
        temp_init_state = np.zeros_like(iter_hmm.init_state) + min_init_val
        temp_init_state[np.clip(np.argmin(np.abs(iter_hmm.gs - mean_freq)), 1, iter_hmm.N - 2)] = 1.
        iter_hmm.init_state = temp_init_state/np.sum(temp_init_state)
    else:
        iter_hmm.init_state = iter_hmm.init_init_state
    return iter_hmm.compute_one_s(loc, tol, max_iter, min_init_val=min_init_val, min_ic=min_ic, save_history=save_history)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path", type=argparse.FileType("rb"), help="path to input dataset")
    parser.add_argument("output_path", type=str, help="base path of output file - appropriate suffixes added later")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--time_before_present", action="store_true", help="dates provided start at a number at the earliest time and count down towards the present")
    group.add_argument("--time_after_zero", action="store_true", help="dates provided start at zero at the earliest time and count up towards the present")
    parser.add_argument("-ytg", "--years_to_gen", type=float, default=1, help="years per generation in VCF or CSV")
    parser.add_argument("-maf", "--min_allele_freq", type=float, default=0.05, help="filters out replicates with mean minor allele frequency < MAF")
    parser.add_argument("--min_sample_density", type=float, default=0.1, help="filters out replicates with fewer than (min_sample_density * max_samples) total samples")
    parser.add_argument("-nc", "--num_cores", type=int, default=1, help="number of CPU cores to parallelize over")
    parser.add_argument("-hs", "--hidden_states", type=int, help="number of approx states in HMM", default=500)
    parser.add_argument("--s_init", type=float, nargs=2, default=[0., 0.], help="vector of initial s value")
    parser.add_argument("-sid", "--starting_init_dist", default="uniform", help="initial initial condition to use")
    parser.add_argument("--sid_dict", nargs='*', help="initial condition dictionary")
    parser.add_argument("-t", "--tol", type=float, default=1e-3, help="ll_i - ll_(i-1) < tol stops the run")
    parser.add_argument("-m", "--maxiter", type=int, default=2000, help="maximum number of iterations")
    parser.add_argument("-Ne", type=int, default=10000, help="effective population size for the HMM")
    parser.add_argument("--hidden_interp", default="chebyshev", help="interpolation of the hidden states (linear vs. Chebyshev nodes for now)")
    parser.add_argument("--ic_update_type", default="beta", help="type of init cond estimation")
    parser.add_argument("--progressbar", action="store_true", help="adds a tqdm progress bar")
    parser.add_argument("--save_history", action="store_true", help="save s and ll history for in-depth analysis (memory intensive)")
    parser.add_argument("--selection_modes", default="all", nargs='*', help="strings of update types to run")
    parser.add_argument("--min_itercount", type=int, default=5, help="minimum number of EM iterations before terminating")
    parser.add_argument("--min_init_val", type=float, default=1e-8, help="minimum value of an init state probability")
    parser.add_argument("--save_csv", action="store_true", help="if inputting a VCF, save a CSV to future reduce pre-processing time")
    parser.add_argument("--info_file", type=argparse.FileType("rb"), help="sample times file (if input = VCF)")
    parser.add_argument("--info_cols", type=str, nargs=2, default=["Genetic_ID","Date_mean"], help="names of the ID and dates columns in the sample times file (if input = VCF)")
    parser.add_argument("--full_output", action="store_true", help="save a pickle file with a full set of outputs (in addition to the CSV)")
    parser.add_argument("--force", type=str, nargs=1, help="if the VCF file only contains homozygous loci, force it to be read as either haploid or diploid")
    parser.add_argument("--no_neutral", action="store_true", help="override the requirement that neutral be run")
    parser.add_argument("--compute_cond", action="store_true", help="compute log P(segregating in sample) for conditionalizing")
    args = parser.parse_args()

    hmm_dd = {}
    hmm_dd["approx_states"] = args.hidden_states
    hmm_dd["s_init"] = args.s_init
    hmm_dd["init_cond"] = args.starting_init_dist
    hmm_dd["tol"] = args.tol
    hmm_dd["ytg"] = args.years_to_gen
    hmm_dd["max_iter"] = args.maxiter
    hmm_dd["Ne"] = args.Ne
    hmm_dd["hidden_interp"] = args.hidden_interp
    hmm_dd["ic_update_type"] = args.ic_update_type
    hmm_dd["selection_modes"] = args.selection_modes
    hmm_dd["sid_dict"] = {}
    hmm_dd["min_ic"] = args.min_itercount
    hmm_dd["min_init_val"] = args.min_init_val
    hmm_dd["maf_thresh"] = args.min_allele_freq
    hmm_dd["min_sample_density"] = args.min_sample_density

    if args.sid_dict is not None:
        for ic_pair in args.sid_dict:
            k, v = ic_pair.split('=')
            try:
                hmm_dd["sid_dict"][k] = float(v)
            except:
                hmm_dd["sid_dict"][k] = v


    hmm_basepath = args.output_path
    pd_path = Path(args.input_path.name)
    num_cores = args.num_cores

    if not pd_path.is_file():
        print(f"no params dict: {pd_path}")
        raise ValueError

    hmm_data = {}
    if pd_path.suffix == ".csv":
        pf = np.loadtxt(pd_path, delimiter="\t")
        hmm_data["final_data"] = pf[:, 2::3].astype(int)
        hmm_data["num_samples"] = pf[:, 1::3].astype(int)
        hmm_data["sample_times"] = (pf[:, ::3] / hmm_dd["ytg"]).astype(int)
        if args.time_before_present:
            hmm_data["sample_times"] = np.max(hmm_data["sample_times"], axis=1)-hmm_data["sample_times"]
        max_samples = np.max(np.sum(hmm_data["num_samples"], axis=1))
    elif pd_path.suffix == ".vcf" or pd_path.suffixes == [".vcf", ".gz"]:
        print("processing VCF...")
        if not (args.info_file and args.info_cols):
            raise ValueError("No sample file (or columns therein) specified!")
        vcf_file = allel.read_vcf(str(pd_path))
        max_samples = vcf_file["calldata/GT"].shape[1] + np.any(vcf_file['calldata/GT'][:, :, 1] >= 0, axis=0).sum()
        vcf_dates = read_csv(args.info_file.name, usecols=args.info_cols, sep="\t").to_numpy()
        if args.time_after_zero:
            vcf_dates[:, 1] = np.max(vcf_dates[:, 1]) - vcf_dates[:, 1]
        full_array = vcf_to_useful_format(vcf_file, vcf_dates, years_per_gen = hmm_dd["ytg"], force=args.force)
        hmm_data["final_data"] = full_array[:, 2::3].astype(int)
        hmm_data["num_samples"] = full_array[:, 1::3].astype(int)
        hmm_data["sample_times"] = full_array[:, ::3].astype(int)
        print("done processing VCF")
    else:
        raise ValueError("File must be a .csv, .vcf, or .vcf.gz!")

    for r_i in np.arange(hmm_data["sample_times"].shape[0]):
        if np.unique(hmm_data["sample_times"][r_i, :]).shape[0] != hmm_data["sample_times"].shape[1]:
            raise ValueError("sampling times must be unique!")

    total_fd = np.sum(hmm_data["final_data"], axis=1)
    total_ns = np.sum(hmm_data["num_samples"], axis=1)

    min_fd = np.minimum(total_fd, total_ns - total_fd)
    MAF_filter_mask = min_fd > total_ns * hmm_dd["maf_thresh"]

    num_samples_mask = np.sum(hmm_data["num_samples"] != 0, axis=1) > 1
    anc_samples_mask = np.sum(hmm_data["num_samples"], axis=1) > hmm_dd["min_sample_density"] * max_samples
    combo_mask = num_samples_mask & MAF_filter_mask & anc_samples_mask

    print(f"MAF filtered {hmm_data['final_data'].shape[0] - MAF_filter_mask.sum()}\n"
          f"< 2 samples {hmm_data['final_data'].shape[0] - num_samples_mask.sum()}\n"
          f"few samples filtered {hmm_data['final_data'].shape[0] - anc_samples_mask.sum()}\n"
          f"overall filtered {hmm_data['final_data'].shape[0] - combo_mask.sum()}")

    hmm_data["final_data"] = hmm_data["final_data"][combo_mask, :]
    hmm_data["num_samples"] = hmm_data["num_samples"][combo_mask, :]
    hmm_data["sample_times"] = hmm_data["sample_times"][combo_mask, :]
    hmm_dd["max_samples"] = max_samples
    hmm_dd["samples_mask"] = combo_mask
    hmm_dd["samples_idxs"] = np.where(combo_mask)[0]

    if args.save_csv:
        csv_path = ""
        if pd_path.suffix == ".vcf":
            csv_path = pd_path.with_suffix(".csv")
        elif pd_path.suffixes == [".vcf", ".gz"]:
            csv_path = pd_path.with_suffix('').with_suffix(".csv")
        if csv_path:
            np.savetxt(csv_path, full_array[combo_mask, :], delimiter="\t", fmt="%d",
                   header="Each row = one replicate; each set of three columns = (sampling time; total samples; derived alleles)")
            print(f"Saved .csv to {csv_path}!")

    if pd_path.suffix == ".vcf" and args.full_output:
        hmm_dd["pos"] = vcf_file["variants/POS"][combo_mask]
        hmm_dd["snp_ids"] = vcf_file["variants/ID"][combo_mask]
        hmm_dd["ref_allele"] = vcf_file["variants/REF"][combo_mask]
        hmm_dd["alt_allele"] = vcf_file["variants/ALT"][combo_mask, 0]

    all_selection_modes = ["neutral", "add", "dom", "rec", "het", "full"]
    if hmm_dd["selection_modes"] == "all" or hmm_dd["selection_modes"] == ["all"]:
        selection_modes = all_selection_modes
        hmm_dd["selection_modes"] = all_selection_modes
    else:
        if len([update_i for update_i in hmm_dd["selection_modes"] if update_i not in all_selection_modes]) > 0:
            raise ValueError("Invalid update type specified!")
        if "neutral" not in hmm_dd["selection_modes"] and not args.no_neutral:
            hmm_dd["selection_modes"] = ["neutral", *hmm_dd["selection_modes"]]
        if "neutral" in hmm_dd["selection_modes"] and hmm_dd["selection_modes"][0] != "neutral":
            hmm_dd["selection_modes"].insert(0,hmm_dd["selection_modes"].pop(hmm_dd["selection_modes"].index("neutral")))
        selection_modes = hmm_dd["selection_modes"]
    for selmode_i, sel_type in enumerate(selection_modes):
        print(f"Now analyzing: {sel_type}!")
        if num_cores > 1:
            parallel_loop = tqdm(range(hmm_data["final_data"].shape[0])) if args.progressbar else range(hmm_data["final_data"].shape[0])
            with Parallel(n_jobs=num_cores) as parallel:
                iter_hmm = HMM(hmm_dd["approx_states"], hmm_dd["Ne"],hmm_dd["s_init"],
                            init_cond = hmm_dd["init_cond"], hidden_interp = hmm_dd["hidden_interp"], **hmm_dd["sid_dict"])
                iter_hmm.update_type = sel_type
                iter_hmm.update_func, iter_hmm.update_func_args = iter_hmm.get_update_func(sel_type, {})
                iter_hmm.init_update_type = hmm_dd["ic_update_type"]
                iter_hmm.init_update_func, iter_hmm.init_params_to_state_func, iter_hmm.init_update_size = iter_hmm.get_init_update_func(hmm_dd["ic_update_type"])
                res = parallel(delayed(run_one_s)(iter_hmm, hmm_data["final_data"][i], hmm_data["num_samples"][i], hmm_data["sample_times"][i], i, hmm_dd["tol"], hmm_dd["max_iter"], hmm_dd["init_cond"] == "data_mean", hmm_dd["min_init_val"], hmm_dd["min_ic"], args.save_history) for i in parallel_loop)

            #silly fix for SLURM issues
            true_rp3 = [rp[3] for rp in res]
            for r_i, rp in enumerate(res):
                if isinstance(rp[3], np.ndarray):
                    if rp[3].shape == (1,):
                        true_rp3[r_i] = rp[3][0]
            hmm_dict = {
                "s_hist": np.array([rp[0] for rp in res]).T.squeeze(),
                "s_final": np.array([rp[1] for rp in res]).T.squeeze(),
                "ll_hist": np.array([rp[2] for rp in res]).T.squeeze(),
                "ic_dist": np.array([true_rp3]).T.squeeze(),
                "itercount_hist": np.array([rp[4] for rp in res]).squeeze(),
                "exit_codes": np.array([rp[5] for rp in res]).squeeze(),
            }
        else:
            iter_hmm = HMM(hmm_dd["approx_states"], hmm_dd["Ne"], np.array(hmm_dd["s_init"]),init_cond=hmm_dd["init_cond"], hidden_interp=hmm_dd["hidden_interp"], **hmm_dd["sid_dict"])
            iter_hmm.update_type = sel_type
            iter_hmm.update_func, iter_hmm.update_func_args = iter_hmm.get_update_func(sel_type, {})
            s_hist, s_final, ll_hist, ic_dist, itercount_hist, exit_codes = iter_hmm.compute_s(hmm_data["final_data"], hmm_data["num_samples"], hmm_data["sample_times"],
                                                                        sel_type, hmm_dd["ic_update_type"], hmm_dd["tol"], hmm_dd["max_iter"], progressbar=args.progressbar, data_mean= hmm_dd["init_cond"] == "data_mean", save_history = args.save_history)
            hmm_dict = {
                "s_hist": s_hist,
                "s_final": s_final,
                "ll_hist": ll_hist,
                "ic_dist": ic_dist,
                "itercount_hist": itercount_hist,
                "exit_codes": exit_codes
            }
        if args.compute_cond and sel_type == "neutral":
            data_matrix = np.zeros((len(hmm_data["final_data"]), len(hmm_data["final_data"][0])*3), dtype=int)
            for i in range(len(hmm_data["final_data"])):
                data_matrix[i, ::3] = hmm_data["sample_times"][i]
                data_matrix[i, 1::3] = hmm_data["num_samples"][i]
                data_matrix[i, 2::3] = hmm_data["final_data"][i]
            zeros_dm = np.copy(data_matrix)
            zeros_dm[:, 2::3] = 0
            ns_dm = np.copy(data_matrix)
            ns_dm[:, 2::3] = ns_dm[:, 1::3]
            if hmm_dd["ic_update_type"] != "fixed":
                est_ic_init_state = np.zeros((data_matrix.shape[0], iter_hmm.gs.shape[0]))
                for i in np.arange(data_matrix.shape[0]):
                    beta_distrib = beta(hmm_dict["ic_dist"][i, 0], hmm_dict["ic_dist"][i, 1])
                    beta_pdf = beta_distrib.pdf(iter_hmm.gs[1:-1])
                    est_ic_init_state[i, 1:-1] = beta_pdf / np.sum(beta_pdf)
            else:
                est_ic_init_state = None
            zeros_lls = iter_hmm.compute_multiple_ll(0,0,zeros_dm,init_states=est_ic_init_state)
            ns_lls = iter_hmm.compute_multiple_ll(0,0,ns_dm,init_states=est_ic_init_state)
            hmm_dict["cond_correction_ll"] = np.log(1-np.exp(zeros_lls)-np.exp(ns_lls))
        
        if args.save_history:
            hmm_dict["ll_final"] = np.array([hmm_dict["ll_hist"][hmm_dict["itercount_hist"][i], i] for i in range(hmm_dict["itercount_hist"].shape[0])])
        else:
            hmm_dict["ll_final"] = hmm_dict["ll_hist"]
            del hmm_dict["ll_hist"]
            del hmm_dict["s_hist"]
        if sel_type == "neutral":
            hmm_dd["neutral_ll"] = hmm_dict["ll_final"]
            hmm_dd["neutral_ic"] = hmm_dict["ic_dist"]
            hmm_dd["neutral_itercount"] = hmm_dict["itercount_hist"]
            if args.compute_cond:
                hmm_dd["cond_correction_ll"] = hmm_dict["cond_correction_ll"]
        else:
            hmm_dd[f"{sel_type}_run"] = hmm_dict

    print("saving output...")
    save_csv_output(hmm_dd, hmm_basepath+".csv")
    if args.full_output:
        with open(hmm_basepath+".pkl", "wb") as file:
            pickle.dump(hmm_dd, file)

if __name__ == "__main__":
    main()
