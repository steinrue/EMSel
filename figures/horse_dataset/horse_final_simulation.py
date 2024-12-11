import matplotlib.pyplot as plt
from pathlib import Path
from cycler import cycler
import numpy as np
from emsel.emsel_util import generate_multiple_wf_data
from tqdm import tqdm
###### MODIFY

data_dir = "data"
output_dir = "output"

###### DO NOT MODIFY

plt.rcParams.update({'font.size': 9,
                     'text.usetex': False,
                     'font.family': 'serif',
                     'font.serif': 'cmr10',
                     'mathtext.fontset': 'cm',
                     'axes.unicode_minus': False,
                     'axes.formatter.use_mathtext': True,})
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
coolormap = plt.get_cmap("Dark2")
colorlist = ["#1D6996", *[coolormap(i) for i in [1,0]], colors[3], colors[4], colors[5]]
plt.rcParams["axes.prop_cycle"] = cycler(color=colorlist)

num_sims = 10000
num_reps = 1000
maf_thresh = .05

for Ne in [16000]:
    for s_i, suffix in enumerate(["full"]):
        for init_est in ["neut", "uncon"]:
            print("-----------")
            seed = 5 + s_i
            sampling_matrix = np.loadtxt(f"{data_dir}/{suffix}_sample_sizes_8.table", skiprows=1, dtype=int)
            sampling_matrix[:, 0] = sampling_matrix[-1, 0] - sampling_matrix[:, 0]
            sampling_matrix = np.flip(sampling_matrix, axis=0)

            nt = sampling_matrix[:, 1]
            sample_times = sampling_matrix.shape[0]
            sample_locs = sampling_matrix[:, 0]
            num_gens = sampling_matrix[-1, 0] + 1

            ps_table = np.loadtxt(f"{data_dir}/{suffix}_{init_est}_means.txt")
            assert ps_table.shape[0] == 2
            p = np.zeros(num_sims) + ps_table[-1]

            final_samples = np.zeros((1, sample_locs.shape[0]))
            final_true_data = np.zeros((1, num_gens))

            i = 0
            while final_samples.shape[0] < num_reps+1:
                temp_samples, temp_true_data = generate_multiple_wf_data(p, 2 * Ne, nt, 0, 0,
                                    num_gens, sample_locs, 5+i, small_s=True)


                # conditioning stuff

                summed_samples = np.sum(temp_samples, axis=1)
                summed_nts = np.sum(nt)
                minor_allele_counts = np.minimum(summed_samples, summed_nts-summed_samples)

                seg_mask = (temp_samples > 0) & (temp_samples < nt)
                summed_seg_times = np.sum(seg_mask, axis=1)
                one_tp_mask = summed_seg_times > 0
                two_tp_mask = summed_seg_times > 1
                maf_mask = minor_allele_counts/summed_nts > maf_thresh

                provisional_samples = temp_samples[two_tp_mask & maf_mask]
                provisional_data = temp_true_data[two_tp_mask & maf_mask]
                print(f"{i}: {provisional_samples.shape[0]}")
                final_samples = np.concatenate((final_samples, provisional_samples))
                final_true_data = np.concatenate((final_true_data, provisional_data))
                i += 1
            final_samples = final_samples[1:num_reps+1, :]
            final_true_data = final_true_data[1:num_reps+1, :]
            print(final_samples.shape)
            data_csv = np.zeros((final_samples.shape[0], nt.shape[0] * 3))
            data_csv[:, ::3] = sample_locs
            data_csv[:, 1::3] = nt
            data_csv[:, 2::3] = final_samples

            data_filename = f"{data_dir}/{suffix}_Ne{Ne}_{init_est}_final_sims.csv"
            np.savetxt(data_filename, data_csv, delimiter="\t", fmt="%d",
                       header="Each row = one replicate; each set of three columns = (sampling time; total samples; derived alleles)")

            fig, axs = plt.subplots(1, 1, figsize=(6, 6), layout="constrained")
            overall_mean = np.median(final_true_data, axis=0)
            upper_q = np.quantile(final_true_data, .95, axis=0)
            lower_q = np.quantile(final_true_data, .05, axis=0)
            # axs.plot(data_subset[:20, :].T, lw=.75)
            axs.plot(overall_mean, color="k", lw=2)
            axs.fill_between(np.arange(overall_mean.shape[0]), lower_q, upper_q, alpha=.5)
            axs.set_title(f"{suffix}_Ne{Ne}_{init_est}")
            axs.set_ylim([0, 1])
            fig.savefig(f"{output_dir}/{suffix}_Ne{Ne}_{init_est}_final_trajs.pdf", format="pdf", bbox_inches="tight")
            plt.close(fig)