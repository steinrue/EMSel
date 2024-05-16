import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import StrMethodFormatter
from numba import njit
from scipy.stats import chi2, beta
from copy import deepcopy

def generate_wf_data(p: float, N: int, N_samples: int, s1: int, s2: int, num_gens: int, sample_locs, seed: int):
    rng = np.random.default_rng(seed=seed)
    freqs = np.zeros(num_gens)
    nt = np.zeros(int(num_gens), dtype=int)
    nt[sample_locs] = N_samples
    freqs[0] = p
    for i in np.arange(num_gens - 1):
        freqs[i + 1] = forward_one_gen(freqs[i], N, s1, s2, rng)
        if freqs[i+1] <= 0:
            all_obs_counts = rng.binomial(nt, freqs)
            return all_obs_counts[sample_locs], freqs
        elif freqs[i+1] >= 1:
            freqs[i+1:] = 1
            all_obs_counts = rng.binomial(nt, freqs)
            return all_obs_counts[sample_locs], freqs

    all_obs_counts = rng.binomial(nt, freqs)
    samples = all_obs_counts[sample_locs]
    true_data = freqs
    return samples, true_data

def generate_multiple_wf_data(p, N, N_samples, s1, s2, num_gens, sample_locs, seed, small_s=False):
    rng = np.random.default_rng(seed=seed)
    freqs = np.zeros((p.shape[0], num_gens))
    freqs[:, 0] = p
    nt = np.zeros_like(sample_locs)+N_samples
    for i in np.arange(num_gens - 1):
        freqs[:, i + 1] = forward_one_gen(freqs[:, i], N, s1, s2, rng, small_s)
        if np.all(freqs[:, i+1] <= 0):
            all_obs_counts = rng.binomial(nt, freqs[:, sample_locs])
            return all_obs_counts, freqs
        elif np.all(freqs[:, i+1] >= 1):
            freqs[:, i+1:] = 1
            all_obs_counts = rng.binomial(nt, freqs[:, sample_locs])
            return all_obs_counts, freqs
    all_obs_counts = rng.binomial(nt, freqs[:, sample_locs])
    return all_obs_counts, freqs

def generate_states_new(n_total, Ne, hidden_interp):
    if hidden_interp == "chebyshev":
        chebyshev_pts = 1/2 + np.cos((2*np.arange(1,n_total-1)-1)*np.pi/(2*(n_total-2)))/2
        all_pts = np.concatenate((np.array([0]), chebyshev_pts[::-1], np.array([1])))
        return all_pts, generate_bounds(all_pts)
    if hidden_interp == "linear":
        return np.linspace(0,1,n_total), generate_bounds(np.linspace(0,1,n_total))
    raise TypeError("Invalid hidden interpolation function!")

def generate_bounds(states):
    bounds = np.zeros(states.shape[0]+1)
    bounds[1:-1] = states[:-1] + np.diff(states)/2
    bounds[0] = states[0]-(states[1]-states[0])/2
    bounds[-1] = states[-1] + (states[-1]-states[-2])/2
    return bounds

def forward_one_gen(p_vector, pop_size, s1, s2, rng, small_s=False):
    if small_s:
        p_prime_vector = p_vector + p_vector * (1 - p_vector) * ((1 - 2 * p_vector) * s1 + p_vector * s2)
    else:
        p_prime_vector = p_vector*(1+s1-s1*p_vector+s2*p_vector)/(1+2*s1*p_vector+s2*p_vector**2-2*s1*p_vector**2)

    return rng.binomial(pop_size, p_prime_vector) / pop_size


def generate_data(pd):
    if "sampling_matrix" in pd:
        nt = pd["sampling_matrix"][:, 1]
        sample_times = pd["sampling_matrix"].shape[0]
        sample_locs = pd["sampling_matrix"][:, 0]
    else:
        nt = np.zeros(int(pd["sample_times"]), dtype=int) + pd["num_samples"]
        sample_times = pd["sample_times"]
        sample_locs = np.linspace(0, int(pd["num_gens"]) - 1, sample_times, dtype=int)
    full_nts = [] if "sampling_matrix" not in pd else np.zeros((1, sample_times))
    full_num_samples = []
    full_sample_locs = []
    p_inits = np.zeros(1,)
    full_samples = np.zeros((1, sample_times))
    full_true_data = np.zeros((1, pd["num_gens"]))
    samples_per_run = 10000
    if pd["init_cond"] == "recip":
        weights = 1/np.arange(1, 2*pd["Ne"])
        weights /= np.sum(weights)
    trial_num = 0
    while full_samples.shape[0]-1 < pd["num_sims"]:
        if pd["init_cond"] == "beta":
            p = np.random.default_rng(pd["seed"]+trial_num).beta(4*pd["Ne"]*pd["mu"], 4*pd["Ne"]*pd["mu"], samples_per_run)
        elif pd["init_cond"] == "fbeta":
            p = np.random.default_rng(pd["seed"]+trial_num).beta(4*pd["Ne"]*pd["mu"], 4*pd["Ne"]*pd["mu"], samples_per_run)
            p[p>.5] = 1-p[p>.5]
        elif pd["init_cond"] == "delta":
            p = np.zeros(samples_per_run)+pd["p_init"]
        elif pd["init_cond"] == "real_special":
            p = np.random.default_rng(pd["seed"]+trial_num).choice(pd["means_array"][1:], size=samples_per_run, replace=True)
        elif pd["init_cond"] == "recip":
            p = np.random.default_rng(pd["seed"]+trial_num).choice(np.arange(1, 2*pd["Ne"])/(2*pd["Ne"]), size=samples_per_run, p=weights)
        else:
            print("weird IC")
            p = pd["p_init"]
        # print(trial_num)
        temp_samples, temp_true_data = generate_multiple_wf_data(p, 2*pd["Ne"], nt, pd["s1_true"], pd["s2_true"],
                                                           pd["num_gens"], sample_locs, pd["seed"] + trial_num, small_s=pd["small_s"])
        if pd["survive_only"]:
            if "sampling_matrix" in pd:
                if "missingness_array" in pd:
                    outer_rng = np.random.default_rng(pd["seed"]+trial_num)
                    missingness_vals = outer_rng.choice(pd["missingness_array"][1:], size=temp_samples.shape[0], replace=True)
                    hits_missing = outer_rng.binomial(temp_samples, 1 - missingness_vals[:, np.newaxis],size=temp_samples.shape)
                    misses_missing = outer_rng.binomial(nt - temp_samples,1 - missingness_vals[:, np.newaxis],size=temp_samples.shape)
                    assert np.all(hits_missing + misses_missing <= nt)
                    print(f"Pre-missingness: {np.mean(temp_samples):.4f} avg. samples. Post: {np.mean(hits_missing):.4f}")

                    temp_nts = hits_missing + misses_missing
                    temp_samples_ms = hits_missing

                    num_samples_mask = np.sum(temp_nts != 0, axis=1) > 1.
                    anc_samples_mask = np.sum(temp_nts, axis=1) > pd["missingness_array"][0] * np.sum(pd["sampling_matrix"][:, 1])
                    total_fd = np.sum(temp_samples_ms, axis=1)
                    total_ns = np.sum(temp_nts, axis=1)

                    min_fd = np.minimum(total_fd, total_ns - total_fd)
                    maf_mask = min_fd > total_ns * pd["means_array"][0]
                    all_mask = anc_samples_mask & num_samples_mask & maf_mask

                    full_nts = np.vstack((full_nts, temp_nts[all_mask, :]))
                    temp_samples = temp_samples_ms
                else:
                    total_fd = np.sum(nt)
                    total_ns = np.sum(temp_samples, axis=1)

                    min_fd = np.minimum(total_fd, total_ns - total_fd)
                    maf_mask = min_fd > total_ns * .05
                    all_mask = maf_mask
            elif pd["sel_type"] == "under":
                samples_greater_than_zero_start = temp_samples[:, 0] > 0
                samples_less_than_one_start = temp_samples[:, 0] < nt[0]
                all_mask = samples_greater_than_zero_start & samples_less_than_one_start
            elif pd["sel_type"] == "over":
                samples_greater_than_zero_endpt = temp_samples[:, -1] > 0
                samples_less_than_all_endpt = temp_samples[:, -1] < nt[-1]
                all_mask = samples_greater_than_zero_endpt & samples_less_than_all_endpt
            else:
                samples_greater_than_zero_endpt = temp_samples[:, -1] > 0
                samples_less_than_all_start = temp_samples[:, 0] < nt[0]
                all_mask = samples_greater_than_zero_endpt & samples_less_than_all_start

            full_samples = np.vstack((full_samples, temp_samples[all_mask, :]))
            full_true_data = np.vstack((full_true_data, temp_true_data[all_mask, :]))
            p_inits = np.hstack((p_inits, p[all_mask]))
        else:
            full_samples = np.vstack((full_samples, temp_samples))
            full_true_data = np.vstack((full_true_data, temp_true_data))
            p_inits = np.hstack((p_inits, p))
        trial_num += 1

        if full_samples.shape[0] < 10 and trial_num * samples_per_run > 1000000:
            break

        if trial_num*samples_per_run > 10000000:
            break

    data_dict = {
        "obs_counts": full_samples[1:pd["num_sims"]+1, :],
        "true_data": full_true_data[1:pd["num_sims"]+1, :],
        "nt": full_nts[1:pd["num_sims"]+1, :] if "sampling_matrix" in pd else nt,
        "sample_times": full_num_samples if full_num_samples else sample_times,
        "sample_locs": full_sample_locs if full_sample_locs else sample_locs,
        "p_inits": p_inits[1:pd["num_sims"]+1],
        "num_runs": trial_num*samples_per_run,
    }

    return data_dict

def get_sel_coeffs_from_type(sel_type, sel_str):
    if sel_type == "add":
        s1 = sel_str / 2
        s2 = sel_str
    elif sel_type == "dom":
        s1 = sel_str
        s2 = sel_str
    elif sel_type == "rec":
        s1 = 0
        s2 = sel_str
    elif sel_type == "over":
        s2 = 0
        s1 = sel_str
    elif sel_type == "under":
        s1 = -sel_str
        s2 = 0
    elif sel_type == "full":
        s1 = sel_str
        s2 = -sel_str
    elif sel_type == "neutral":
        s1 = 0
        s2 = 0
    elif sel_type == "weird":
        s1 = sel_str
        s2 = -2*sel_str
    else:
        raise ValueError("Invalid selection type!")
    return s1, s2

def plot_one_ll_grid(fig, axs, s1_grid, s2_grid, ll_grid, vmin, vmax, s1_label = False, s2_label = False, plot_max = True):
    if vmin == vmax:
        vmin = np.amin(ll_grid)
        vmax = np.amax(ll_grid)
    gridplot_norm = axs.contourf(s1_grid, s2_grid, ll_grid.T,levels=25, cmap="viridis", vmin=vmin, vmax=vmax)
    fig.colorbar(gridplot_norm, ax=axs, shrink=.3, label="Log likelihood")
    axs.axvline(c="k", ls="--", linewidth=.7)
    axs.axhline(c="k", ls="--", linewidth=.7)
    axs.set_xlim(s1_grid[0], s1_grid[-1])
    axs.set_ylim(s2_grid[0], s2_grid[-1])
    if plot_max:
        s1_max_idx_norm, s2_max_idx_norm = np.unravel_index(np.argmax(ll_grid, axis=None), ll_grid.shape)
        axs.plot(s1_grid[s1_max_idx_norm], s2_grid[s2_max_idx_norm], "*",color="blue", ms=6, label="Max.", markeredgewidth=.5, markeredgecolor="k")
    if s1_label:
        axs.set_xlabel("$s_1$")
    if s2_label:
        axs.set_ylabel("$s_2$")
    axs.set_aspect("equal")

def plot_qq(axs, axins, logps, labels, colors=None, legend_loc="upper right", thin=False):
    assert len(logps) == len(labels)
    len_ps = logps[0].shape[0]
    xrange = np.arange(1, len_ps+1)
    num_conf_pts = 500

    max_y = -np.log10(1/len_ps)
    for i in range(len(logps)):
        assert logps[i].shape[0] == len_ps

    for i in range(len(logps)):
        sorted_logps = np.sort(logps[i])[::-1]
        max_y = max(max_y, np.max(logps[i]))
        if thin:
            rng = np.random.default_rng(i).uniform(size=len_ps)
            idxs = np.where((rng > .9) | (sorted_logps > 1))[0]
        else:
            idxs = np.arange(len_ps)
        if colors:
            axins.plot(xrange[idxs] / len_ps, np.power(10, -sorted_logps)[idxs], lw=2, label=labels[i], color=colors[i])
            axs.plot(-np.log10(xrange[idxs] / len_ps), sorted_logps[idxs], ".", lw=2, color=colors[i])
        else:
            axins.plot(xrange[idxs] / len_ps, np.power(10, -sorted_logps)[idxs], lw=2, label=labels[i])
            axs.plot(-np.log10(xrange[idxs] / len_ps), sorted_logps[idxs], ".", lw=2)

    unlog_max_y = np.power(10, -max_y)
    lin_xspace = np.linspace(unlog_max_y/2, len_ps+.9, num_conf_pts)
    geom_xspace = np.geomspace(unlog_max_y/2, len_ps+.9, num_conf_pts)
    lin_sigma = np.sqrt((lin_xspace*(len_ps-lin_xspace+1))/((len_ps+2)*(len_ps+1)**2))
    geom_sigma = np.sqrt((geom_xspace*(len_ps-geom_xspace+1))/((len_ps+2)*(len_ps+1)**2))

    axins.plot([0, 1], [0, 1], color="k", lw=1)
    axins.fill_between(lin_xspace/len_ps, lin_xspace/len_ps-1.96*lin_sigma, lin_xspace/len_ps+1.96*lin_sigma, color="k", alpha=.2)
    axins.set_xlim([0, 1])
    axins.set_ylim([0, 1])
    axs.plot([0, max_y * 1.05], [0, max_y * 1.05], color="k", lw=1)
    upper_loglim = -np.log10(np.array([beta.ppf(.975, i, len_ps-i) for i in range(1,len_ps)]))
    lower_loglim = -np.log10(np.array([beta.ppf(.025, i, len_ps-i) for i in range(1,len_ps)]))
    conf_xspace = -np.log10((np.arange(1,len_ps)-.5)/len_ps)
    #lower_loglim = -np.log10(np.clip(geom_xspace/len_ps+1.96*geom_sigma, np.power(10, -max_y*2), None))
    axs.fill_between(conf_xspace, lower_loglim, upper_loglim, color="k", alpha=.08)


    axs.get_xaxis().set_major_formatter(StrMethodFormatter('{x:,.1f}'))
    axs.get_yaxis().set_major_formatter(StrMethodFormatter('{x:,.1f}'))
    axins.get_xaxis().set_major_formatter(StrMethodFormatter('{x:,.1f}'))
    axins.get_yaxis().set_major_formatter(StrMethodFormatter('{x:,.1f}'))
    #axs.set_yticks(np.arange(max_y*1.05))
    #a#xs.set_yticks(axs.get_yticks()[1:])
    #axs.get_xticklabels()[0].set_horizontalalignment("right")
    axins.set_xticks([0, 0.5, 1])
    axins.set_yticks([0, 0.5, 1])
    axins.set_yticks(axins.get_yticks()[1:])
    axins.get_xticklabels()[0].set_horizontalalignment("right")

    handles, labels = axins.get_legend_handles_labels()
    axs.legend(handles, labels, loc=legend_loc)
    axs.set_xlabel("$\mathbb{E}(-\log_{10}(p))$")
    axs.set_ylabel("$-\log_{10}(p)$")
    axs.set_ylim([0, max_y * 1.05])
    axs.set_xlim([0, max_y * 1.05])

@njit
def get_uq_a_exps(a, powers):
    a_exps_list = [np.linalg.matrix_power(a, powers[0])]
    powers_diff = np.diff(powers)
    for power_diff in powers_diff:
        if power_diff == 1:
            a_exps_list.append(a_exps_list[-1] @ a)
        else:
            a_exps_list.append(a_exps_list[-1] @ np.linalg.matrix_power(a, power_diff))
    return a_exps_list

def classify_full_run(s_table):
    classified_array = np.zeros(s_table.shape[1])
    for idx in np.arange(s_table.shape[1]):
        if np.all(np.isclose(s_table[:, idx], 0)):
            classified_array[idx] = 3
        elif s_table[0, idx] >= 0 and s_table[1, idx] >= 0 and s_table[1, idx] >= s_table[0, idx]:
            classified_array[idx] = 2
        elif s_table[0, idx] <= 0 and s_table[1, idx] <= 0 and s_table[1, idx] <= s_table[0, idx]:
            classified_array[idx] = 2
        elif s_table[0, idx] >= 0:
            classified_array[idx] = 0
        else:
            classified_array[idx] = 1

    summed_array = [np.sum(classified_array == i) for i in np.arange(4)]
    return classified_array, summed_array

def vcf_to_useful_format(vcf_file, sample_times_file, years_per_gen=28.1, force=None):
    sample_times_ordered = np.copy(sample_times_file)
    sample_times_ordered[:, 1] //= years_per_gen
    max_sample_time = np.max(sample_times_ordered[:, 1])
    sample_times_ordered = sample_times_ordered[np.argsort(sample_times_ordered[:, 0]), :]
    correct_order_idxs = vcf_file["samples"].argsort().argsort()
    sample_times_ordered = sample_times_ordered[correct_order_idxs, :]
    sample_times, sample_idxs = np.unique(sample_times_ordered[:, 1], return_inverse=True)
    chroms = vcf_file["variants/CHROM"].astype(int)
    locus_list_list = []

    #if we're doing genome-wide thresholds?
    diploid_tf = np.any(vcf_file["calldata/GT"][:, :, 1] >= 0, axis=0)
    if np.all(vcf_file["calldata/GT"][:, :, 0] == vcf_file["calldata/GT"][:, :, 1]):
        if not force or force not in ["haploid", "diploid"]:
            raise TypeError("VCF call data is all homozygotes - must use --force [haploid/diploid]!")
        elif force == "haploid":
            vcf_file["calldata/GT"][:, :, 1] = -1
    big_final_table = np.zeros((1, sample_times.shape[0]*3))
    for chrom in np.unique(chroms):
        print(chrom)
        final_table = np.zeros(((chroms == chrom).sum(), sample_times.shape[0]*3))
        for sample_i, sample_time in enumerate(sample_times):
            sample_indices = np.where(sample_i == sample_idxs)[0]
            sample_time_indices = np.where(sample_time == sample_times)[0]
            assert (vcf_file['samples'][sample_indices] == sample_times_ordered[sample_indices, 0]).all()
            relevant_calls_nd = np.squeeze(vcf_file["calldata/GT"][:, sample_indices, :])
            num_samples_nd = np.sum(relevant_calls_nd >= 0, axis=-1)
            num_zeros_nd = np.sum(relevant_calls_nd == 0, axis=-1)
            if relevant_calls_nd.ndim > 2:
                num_samples_nd = np.sum(num_samples_nd, axis=1)
                num_zeros_nd = np.sum(num_zeros_nd, axis=1)
            final_data_nd = num_samples_nd - num_zeros_nd
            final_table[:, sample_i*3+1] = num_samples_nd.astype(int)
            final_table[:, sample_i*3+2] = final_data_nd.astype(int)
        final_table[:, ::3] = (max_sample_time - sample_times[::-1]).astype(int)
        final_table[:, 1::3] = final_table[:, 1::3][:, ::-1]
        final_table[:, 2::3] = final_table[:, 2::3][:, ::-1]
        big_final_table = np.vstack((big_final_table, final_table))
    return big_final_table[1:, :]


def full_bh_procedure(llgka_list, fitted_dist, lr_shift, alpha, bh=True):
    p_vals = []
    for llgka in llgka_list:
        temp_lrs = get_lr_statistic(llgka)
        temp_p_vals = np.zeros_like(temp_lrs)
        if lr_shift > 0:
            temp_p_vals[temp_lrs>lr_shift] = (1-fitted_dist.cdf(temp_lrs[temp_lrs>lr_shift]-lr_shift))/2
            temp_p_vals[temp_lrs<=lr_shift] = np.clip(1-temp_lrs[temp_lrs<=lr_shift]/(2*lr_shift), .5, 1)
        else:
            temp_p_vals = 1-fitted_dist.cdf(temp_lrs)
        p_vals.append(temp_p_vals)

    flat_p_vals = np.zeros(1)
    for p_val_array in p_vals:
        for col in np.arange(p_val_array.shape[1]):
            flat_p_vals = np.hstack((flat_p_vals, p_val_array[:, col]))
    flat_p_vals = flat_p_vals[1:]
    if bh:
        BH_line, rejected_ps = bh_correct(flat_p_vals, alpha, yekutieli=False)
    else:
        BH_line, rejected_ps = (alpha, np.where(flat_p_vals <= alpha)[0])
    p_idx = 0
    classified_array_list = []
    for p_i, array_p in enumerate(p_vals):
        classified_array = np.zeros_like(array_p)
        for col_p in np.arange(classified_array.shape[1]):
            non_neutral_idxs = np.intersect1d(rejected_ps[rejected_ps < p_idx + classified_array.shape[0]], rejected_ps[rejected_ps > p_idx])-p_idx
            if non_neutral_idxs.shape[0] > 0:
                classified_array[non_neutral_idxs, col_p] = np.argmax(llgka_list[p_i][non_neutral_idxs, 1:, col_p], axis=1)+1
            p_idx += classified_array.shape[0]
        classified_array_list.append(classified_array)
    return BH_line, p_vals, classified_array_list

def bh_correct(p_values, alpha, yekutieli=False):
    M = p_values.shape[0]
    p_values_sorted = np.sort(p_values.copy())
    bh_range = np.arange(1,M+1)
    if yekutieli:
        alpha /= np.sum(1/bh_range)
    small_ps = np.where(p_values_sorted <= bh_range*alpha/M)[0]
    if small_ps.shape[0] > 0:
        k_max = np.where(p_values_sorted <= bh_range*alpha/M)[0][-1]
    else:
        print("no significant ps here!")
        return 1, np.array([])
    p_k = np.sqrt(p_values_sorted[k_max]* p_values_sorted[k_max+1])
    return p_k, np.where(p_values <= p_k)[0]

def get_llg_array(one_data_dict, onep_types, full_classified_array):
    ll_array = np.zeros((one_data_dict["neutral_ll"].shape[0], len(onep_types) + 4))
    ll_array[:, 0] = one_data_dict["neutral_ll"]
    for i, onep_type in enumerate(onep_types):
        ll_array[:, i+1] = one_data_dict[f"{onep_type}_run"]["ll_final"]
    ll_array[:, -3] = np.where(full_classified_array == 0, one_data_dict["full_run"]["ll_final"], -np.inf)
    ll_array[:, -2] = np.where(full_classified_array == 1, one_data_dict["full_run"]["ll_final"], -np.inf)
    ll_array[:, -1] = np.where(full_classified_array == 2, one_data_dict["full_run"]["ll_final"], -np.inf)
    return ll_array

def get_llgka_array(llg_array, k, alpha):
    llgka_array = np.copy(llg_array)
    llgka_array[:, 0, :] += alpha
    llgka_array[:, -3:, :] -= k
    return llgka_array

def get_one_count_matrix(k, alt_array):
    altk_array = get_llgka_array(alt_array, k, 0)
    count_matrix = np.zeros((alt_array.shape[2], alt_array.shape[1] - 1))
    for j in np.arange(alt_array.shape[2]):
        count_matrix[j, :] = classify_llgk_array(altk_array[:, 1:, j])
    return count_matrix

def classify_llgk_array(llgk_array):
    max_locs = np.argmax(llgk_array, axis=1)
    min_counts = [np.sum(max_locs == i) for i in np.arange(llgk_array.shape[1])]
    return min_counts

def get_lr_statistic(llgk_array):
    best_lr = np.amax(llgk_array[:, 1:, :], axis=1)
    return 2*(best_lr - llgk_array[:, 0, :])

def params_dict_to_str(**kwargs):
    name_str = ""
    NEUTRAL_FLAG = False
    if "sel_type" in kwargs.keys():
        name_str += kwargs["sel_type"] + "_"
        NEUTRAL_FLAG = kwargs["sel_type"] == "neutral"
    if "sel_str" in kwargs.keys() and not NEUTRAL_FLAG:
        name_str += "s"+ str(kwargs["sel_str"])[2:] + "_"
    if "num_gens" in kwargs.keys():
        name_str += "g" + str(kwargs["num_gens"]) + "_"
    if "init_dist" in kwargs.keys():
        if kwargs["init_dist"] == "fbeta":
            name_str += "fbeta" + "_"
        else:
            name_str += "d" + str(kwargs["init_dist"])[2:] + "_"
    if "condo" in kwargs.keys():
        if kwargs["condo"]:
            name_str += "cond" + "_"
        else:
            name_str += "uncond" + "_"
    if "num_replicates" in kwargs.keys():
        name_str += "nr" + str(kwargs["num_replicates"]) + "_"
    if "sample_times" in kwargs.keys():
        name_str += "st" + str(kwargs["sample_times"]) + "_"
    if "num_samples" in kwargs.keys():
        name_str += "ns" + str(kwargs["num_samples"]) + "_"
    if "Ne" in kwargs.keys():
        name_str += "Ne" + str(kwargs["Ne"]) + "_"
    if "hidden_states" in kwargs.keys():
        name_str += "hs" + str(kwargs["hidden_states"]) + "_"
    return name_str[:-1]

def generate_real_ic(inv_cdf, num_replicates, seed, betas=False):
    rng = np.random.default_rng(seed)
    if betas:
        init_vals = np.zeros(num_replicates)
        beta_choices = rng.choice(inv_cdf.shape[0], size=num_replicates, p=inv_cdf[:, 0])
        for beta_i, beta in enumerate(beta_choices):
            init_vals[beta_i] = rng.beta(inv_cdf[beta, 1], inv_cdf[beta, 2])
        return init_vals
    else:
        return inv_cdf(rng.uniform(size=num_replicates))

def average_p_vals(log_p_vals, positions, window=10000, ev=True):
    avg_p_vals = np.zeros_like(log_p_vals)
    if window > 1000:
        pos_windows = np.zeros((positions.shape[0], 2))
        pos_windows[:, 0] = np.clip(positions-window, a_min=0, a_max=None)
        pos_windows[:, 1] = np.clip(positions+window, a_min=None, a_max=positions[-1])
        pos_locs = np.searchsorted(positions, pos_windows)
        expected_snps = np.mean(np.diff(pos_locs, axis=1))
        for i in np.arange(log_p_vals.shape[0]):
            avg_p_vals[i] = np.sum(log_p_vals[np.arange(pos_locs[i, 0], pos_locs[i, 1])])
        if ev:
            avg_p_vals /= expected_snps
    else:
        for i in np.arange(avg_p_vals.shape[0]):
            avg_p_vals[i] = np.mean(log_p_vals[max(i-window, 0):min(i+window,avg_p_vals.shape[0]-1)])
    return avg_p_vals

def get_roc_and_auc(n_p_vals, nn_p_vals):
    n_denom = n_p_vals.shape[0]
    nn_denom = nn_p_vals.shape[0]
    num_pts = 1000
    TPR = np.zeros(num_pts+2)
    FPR = np.zeros(num_pts+2)
    TPR[0] = FPR[0] = 0
    TPR[-1] = FPR[-1] = 1
    for t_i, thresh in enumerate(np.logspace(-10, 0, num_pts)):

        TPR[t_i+1] = (nn_p_vals < thresh).sum()/nn_denom
        FPR[t_i+1] = (n_p_vals < thresh).sum()/n_denom

    auc = 0
    for i in np.arange(1,num_pts+1):
        auc += 1/2*(TPR[i]-TPR[i-1])*(FPR[i]-FPR[i-1])
        auc += TPR[i]*(FPR[i+1]-FPR[i])

    return FPR, TPR, auc

def constrain_step(s1_0, s2_0, s1_prime, s2_prime, step_constraint):
    if np.sqrt((s1_0 - s1_prime) ** 2 + (s2_0 - s2_prime) ** 2) > step_constraint:
        print("step too large!")
        s_dist = np.sqrt((s1_0 - s1_prime) ** 2 + (s2_0 - s2_prime) ** 2)
        s1_prime = s1_0 + (s1_prime - s1_0) * step_constraint / s_dist
        s2_prime = s2_0 + (s2_prime - s2_0) * step_constraint / s_dist
    return s1_prime, s2_prime

def get_sig_pts(p_vals, prim_logp, sec_logp, window, num_sig):
    p_second_mask = p_vals > sec_logp

    sig_locs = np.where(p_vals > prim_logp)[0]
    sig_pass = []
    for sig_loc in sig_locs:
        num_sig_temp = np.sum(p_second_mask[sig_loc - window:sig_loc + window + 1]) - 1
        if num_sig_temp >= num_sig:
            sig_pass.append(sig_loc)

    sig_idxs = np.array(sig_pass)
    return sig_idxs

def get_1d_s_data_from_type(s_data, sel_type):
    if sel_type in ["add", "rec"]:
        return s_data[1, :]
    elif sel_type in ["dom", "over", "under", "het", "full"]:
        return s_data[0, :]
    else:
        raise ValueError


def convert_from_abbrevs(names_list, shorthet=False, shortall=False):
    short_names = ["add", "dom", "rec", "over", "under", "het", "full"]
    long_names = ["Additive", "Dominant", "Recessive", "Overdominant", "Underdominant", "Heterozygote difference", "Multi-alternative"]
    if shorthet:
        long_names[long_names.index("Heterozygote difference")] = "Het. diff."
    if shortall:
        long_names = ["Add.", "Dom.", "Rec.", "Over.", "Under.", "Het. diff.", "Uncons."]
    if isinstance(names_list, list):
        names_list_copy = deepcopy(names_list)
        for l_i, long_name in enumerate(long_names):
            names_list_copy = [long_name if name_val==short_names[l_i] else name_val for name_val in names_list_copy]
        return names_list_copy
    elif isinstance(names_list, str):
        if names_list in short_names:
            return long_names[short_names.index(names_list)]
        return names_list
    else:
        raise ValueError

def convert_to_abbrevs(names_list):
    short_names = ["add", "dom", "rec", "over", "under", "het", "full"]
    long_names = ["Additive", "Dominant", "Recessive", "Overdominant", "Underdominant", "Heterozygote difference", "Multi-alternative"]
    if isinstance(names_list, list):
        names_list_copy = deepcopy(names_list)
        for s_i, short_name in enumerate(short_names):
            names_list_copy = [short_name if name_val==long_names[s_i] else name_val for name_val in names_list_copy]
        return names_list_copy
    elif isinstance(names_list, str):
        if names_list in long_names:
            return short_names[long_names.index(names_list)]
        return names_list
    else:
        raise ValueError


def extendedFisher(pMatrix, standardBrown=False, fancyBrown=False, estimateLoc=False, epsilon=1e-16):
    (L, n) = pMatrix.shape

    # statistic is 2 times negtive sum of logarithms
    fisherStat = 2 * np.sum(pMatrix, axis=1)

    # browns methods attempt to correct for correlation by taking mean and var of fisherStat into account
    if (standardBrown):
        # get stats of fisherStat
        fisherMean = np.mean(fisherStat)
        fisherVar = np.var(fisherStat)

        # set chi2 parameters according to brown
        theLoc = 0
        theScale = fisherVar / (2 * fisherMean)
        df = (2 * fisherMean * fisherMean) / fisherVar
    elif (fancyBrown):
        # fit parameters
        if (estimateLoc):
            (df, theLoc, theScale) = chi2.fit(fisherStat)
        else:
            (df, theLoc, theScale) = chi2.fit(fisherStat, floc=0)
            # (df, lam, theLoc, theScale) = scipy.stats.ncx2.fit (fisherStat, floc=0)
            # return numpy.clip (1 - scipy.stats.ncx2.cdf(fisherStat, df, lam, scale=theScale, loc=theLoc), epsilon, 1-epsilon)
    else:
        # just regular fisher
        theLoc = 0
        theScale = 1
        df = 2 * n

    # use the respective chi2 for the p-values
    fit_chi2 = chi2(df, scale=theScale, loc=theLoc)
    all_log_p = -fit_chi2.logsf(fisherStat)/np.log(10)

    return all_log_p


def windowMatrix (values, halfWindowSize=25):
    # offset to only those where we have full data
    fullWindowSize = 1 + 2 * halfWindowSize
    extendedData = np.zeros ((len(values)-fullWindowSize, fullWindowSize))
    # go through indices for beginning
    for bIdx in np.arange(extendedData.shape[1]):
        # get corresponding index for end
        eIdx = bIdx + extendedData.shape[0]
        # and copy
        extendedData[:,bIdx] = values[bIdx:eIdx]
    return extendedData

def save_csv_output(hmm_dict, hmm_path):
    num_reps = hmm_dict["neutral_ll"].shape[0]
    info_array = np.zeros((num_reps, 2+3*(len(hmm_dict["selection_modes"])-1)), dtype=float)
    info_array[:, 0] = hmm_dict["samples_idxs"]
    info_array[:, 1] = hmm_dict["neutral_ll"].flatten()
    cols = ["unfiltered index", "neutral_ll"]
    for ud_i, ud_type in enumerate(hmm_dict["selection_modes"][1:], start=1):
        info_array[:, 3*ud_i-1] = hmm_dict[f"{ud_type}_run"]["ll_final"].flatten()
        info_array[:, 3*ud_i:3*(ud_i+1)-1] = hmm_dict[f"{ud_type}_run"]["s_final"].T
        cols.extend([f"{ud_type}_ll", f"{ud_type}_s1", f"{ud_type}_s2"])
    info_df = pd.DataFrame(info_array, columns=cols)
    info_df.to_csv(hmm_path, index=False)