import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from pathlib import Path
import pickle
from emsel.emsel_util import classify_full_run, correct_for_het, get_one_count_matrix, params_dict_to_str, get_llg_array, get_lr_statistic, get_llgka_array, full_bh_procedure, convert_from_abbrevs, convert_to_abbrevs, plot_qq
from copy import deepcopy
from scipy.stats import chi2, gengamma
from pandas import DataFrame
from scipy.signal import savgol_filter
import seaborn as sns
from cycler import cycler
from tqdm import tqdm


###### MODIFY
sel_strs = [.005, .01, .025, .05]
num_gens_list = [101,251,1001]
init_dists = [.005, .25, "recip"]

num_gens_list = [125]
init_dists = ["real_special"]

EM_dir = "EM/ibdne/big"
data_dir = "data/ibdne/big"
output_dir = "output/ibdne"
classified_dir = "classified"

###### DO NOT MODIFY

if "matched" in EM_dir:
    EM_suffix = "Ne10496_"
    output_suffix = "real_matched_"
elif "ibdne" in EM_dir:
    EM_suffix = "Ne35313_"
    output_suffix = "ibdne_"
else:
    EM_suffix = ""
    output_suffix = ""

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
init_colorlist = colorlist
plt.rcParams["axes.prop_cycle"] = cycler(color=colorlist)

classification_types = ["neutral", "add", "dom", "rec", "over", "under", "full"]
classification_types_rows = ["Neut.", "Add.", "Dom.", "Rec.", "Over.", "Under.", "Indet."]
sel_types = ["add", "dom", "rec", "over", "under"]
sel_types_rows = ["add", "dom", "rec", "over", "under"]

onep_types = ["add", "dom", "rec", "het"]
update_types = [*onep_types, "full"]
final_sel_types = ["add", "dom", "rec", "over", "under"]

p_val = .05
k_opt = 1
save_bh = True

neutral_ll_array = np.array([0])
full_ll_array = np.array([0])
uncon_ll_array = np.array([0])
onep_ll_array = np.array([0]*len(onep_types))

bigtable_idxs = []
full_bigtable_list = []
for num_gens in num_gens_list:
    for init_dist in init_dists:
        row_list = []
        g_list_data = []
        llrs_data = []
        hs_data = []
        ndict = {}
        ndict["sel_type"] = "neutral"
        ndict["num_gens"] = num_gens
        ndict["init_dist"] = init_dist

        neutral_filename = params_dict_to_str(**ndict)
        row_list.append([neutral_filename])

        for e_i, suff in enumerate(tqdm(["plot1", "plot2", "plot3"])):
            neutral_EM_path = Path(f"{EM_dir}/{neutral_filename}_{suff}_{EM_suffix}EM.pkl")
            with open(neutral_EM_path, "rb") as file:
                nf_data = pickle.load(file)
            neutral_ll_array = np.concatenate((neutral_ll_array, nf_data["neutral_ll"]))
            full_ll_array = np.concatenate((full_ll_array, nf_data["full_run"]["ll_final"]))
            uncon_ll_array = np.concatenate((uncon_ll_array, nf_data["full_run"]["ll_final"]))
            neutral_llg = get_llg_array(nf_data, onep_types, np.zeros_like(nf_data["neutral_ll"]))
            onep_ll_array = np.vstack((onep_ll_array, neutral_llg[:, 1:len(onep_types)+1]))


            n_data = np.loadtxt(Path(f"{data_dir}/{neutral_filename}_{suff}_data.csv", sep="\t", fmt="%d"))
            if e_i == 0:
                data = n_data
            else:
                data = np.vstack((data, n_data))

        neutral_ll_array = neutral_ll_array[1:]
        uncon_ll_array = uncon_ll_array[1:]
        full_ll_array = full_ll_array[1:]
        full_ll_array -= 1
        onep_ll_array = onep_ll_array[1:]
        full_ll_array = full_ll_array[:, np.newaxis]

        full_ll_array = np.hstack((onep_ll_array, full_ll_array))

        # neutral_data_path = Path(f"{EM_dir}/{neutral_filename}_{EM_suffix}EM.pkl")
        # with open(neutral_data_path, "rb") as file:
        #     nf_data = pickle.load(file)
        #
        #
        # neutral_ll_array = nf_data["neutral_ll"]
        # llg_array = get_llg_array(nf_data, onep_types, np.zeros_like(nf_data["neutral_ll"]))
        # print(llg_array[:10])
        # llg_array[:, len(onep_types)+1:] -= 1
        #
        # onep_ll_array = llg_array[:, 1:len(onep_types)+1]
        # full_ll_array = llg_array[:, 1:]
        #
        # print(full_ll_array.shape)
        # print((np.argmax(full_ll_array, axis=1)==4).sum())
        # print((np.argmax(full_ll_array, axis=1)==5).sum())
        # print((np.argmax(full_ll_array, axis=1)==6).sum())
        # print(num_gens)
        # print(init_dist)
        # raise Error


        finaluncon_gengamma_path = Path(f"{output_dir}/{neutral_filename}_{output_suffix}finaluncon_gengamma_fit.pkl")
        final_gengamma_path = Path(f"{output_dir}/{neutral_filename}_{output_suffix}final_gengamma_fit.pkl")

        finaluncon_chi2_path = Path(f"{output_dir}/{neutral_filename}_{output_suffix}finaluncon_chi2_fit.pkl")
        final_chi2_path = Path(f"{output_dir}/{neutral_filename}_{output_suffix}final_chi2_fit.pkl")

        print(final_gengamma_path)
        gengamma_path = Path(f"{output_dir}/{neutral_filename}_{output_suffix}gengamma_resim_fit.pkl")
        chi2_path = Path(f"{output_dir}/{neutral_filename}_{output_suffix}chi2_resim_fit.pkl")

        with open(finaluncon_gengamma_path, "rb") as file:
            finaluncon_gengamma_dict = pickle.load(file)

        with open(final_gengamma_path, "rb") as file:
            final_gengamma_dict = pickle.load(file)

        with open(finaluncon_chi2_path, "rb") as file:
            finaluncon_chi2_dict = pickle.load(file)

        with open(final_chi2_path, "rb") as file:
            final_chi2_dict = pickle.load(file)

        with open(gengamma_path, "rb") as file:
            gengamma_dict = pickle.load(file)

        with open(chi2_path, "rb") as file:
            chi2_dict = pickle.load(file)

        n_llr_data = 2*(uncon_ll_array-neutral_ll_array)


        # a = np.where(n_llr_data == np.max(n_llr_data))[0]
        # print(a)
        # for val in a:
        #     print(data[val, 2::3])
        #     print(np.sum(data[val]))
        # plt.plot(data[val, ::3], data[val, 2::3]/data[val, 1::3])
        # plt.savefig(Path(f"{output_dir}/{neutral_filename}_reviewefinal_bad_data.pdf"))
        # raise Error
        # med_p_vals = np.zeros_like(n_llr_data)
        # med_p_vals[n_llr_data > gengamma_dict["lr_shift"]] = (1 - gengamma_dict["dist_fit"].cdf(n_llr_data[n_llr_data > gengamma_dict["lr_shift"]] - gengamma_dict["lr_shift"])) / 2
        # med_p_vals[n_llr_data <= gengamma_dict["lr_shift"]] = np.clip(1 - n_llr_data[n_llr_data <= gengamma_dict["lr_shift"]] / (2 * gengamma_dict["lr_shift"]), .5, 1)
        #
        # fit_chi2_p_vals = np.zeros_like(n_llr_data)
        # fit_chi2_p_vals[n_llr_data > chi2_dict["lr_shift"]] = (1 - chi2_dict["dist_fit"].cdf(n_llr_data[n_llr_data > chi2_dict["lr_shift"]] - chi2_dict["lr_shift"])) / 2
        # fit_chi2_p_vals[n_llr_data <= chi2_dict["lr_shift"]] = np.clip(1 - n_llr_data[n_llr_data <= chi2_dict["lr_shift"]] / (2 * chi2_dict["lr_shift"]), .5, 1)
        #
        # p_vals_chi21 = 1 - chi2(1).cdf(n_llr_data)
        # p_vals_chi22 = 1 - chi2(2).cdf(n_llr_data)
        #
        # fig, axs = plt.subplots(1, 1, figsize=(3.1, 3.1), layout="constrained")
        # logps = [-np.log10(med_p_vals), -np.log10(p_vals_chi21), -np.log10(p_vals_chi22),
        #          -np.log10(fit_chi2_p_vals)]
        # labels = ["Gen. gamma", r"$\chi^2(1)$", r"$\chi^2(2)$", rf"$\chi^2(k={{{chi2_dict['dist_fit'].args[0]:.2f}}})$"]
        # colors = [colorlist[2], "#861A5E", colorlist[0], "r"]
        # axins = axs.inset_axes([.65, .11, .3, .3])
        # #axs.text(-.2, .97, r"$\bf{D}$", fontsize=13, transform=axs.transAxes)
        # plot_qq(axs, axins, logps, labels, colors=colors, legend_loc="upper left", thin=True, rasterized=True)
        # fig.savefig(
        #     Path(f"{output_dir}/neutral_g{num_gens}_d{init_dist}_{output_suffix}uncon_100k_d2_newsample.pdf"),
        #     format="pdf", bbox_inches="tight", dpi=300)
        # plt.close(fig)

        llr_data = onep_ll_array
        fit_gengamma = final_gengamma_dict
        fit_chi2 = final_chi2_dict
        title = "final_1dmax"
        final_llr_data = 2*(np.max(llr_data, axis=1)-neutral_ll_array)

        final_med_p_vals = np.zeros_like(final_llr_data)
        final_med_p_vals[final_llr_data > fit_gengamma["lr_shift"]] = (1 - fit_gengamma["dist_fit"].cdf(
            final_llr_data[final_llr_data > fit_gengamma["lr_shift"]] - fit_gengamma["lr_shift"])) / 2
        final_med_p_vals[final_llr_data <= fit_gengamma["lr_shift"]] = np.clip(
            1 - final_llr_data[final_llr_data <= fit_gengamma["lr_shift"]] / (2 * fit_gengamma["lr_shift"]), .5, 1)

        final_fit_chi2_p_vals = np.zeros_like(final_llr_data)
        final_fit_chi2_p_vals[final_llr_data > fit_chi2["lr_shift"]] = (1 - fit_chi2["dist_fit"].cdf(
            final_llr_data[final_llr_data > fit_chi2["lr_shift"]] - fit_chi2["lr_shift"])) / 2
        final_fit_chi2_p_vals[final_llr_data <= fit_chi2["lr_shift"]] = np.clip(
            1 - final_llr_data[final_llr_data <= fit_chi2["lr_shift"]] / (2 * fit_chi2["lr_shift"]), .5, 1)

        final_p_vals_chi21 = 1 - chi2(1).cdf(final_llr_data)
        final_p_vals_chi22 = 1 - chi2(2).cdf(final_llr_data)

        fig, axs = plt.subplots(1, 1, figsize=(3.1, 3.1), layout="constrained")
        logps = [-np.log10(final_p_vals_chi21), -np.log10(final_p_vals_chi22), -np.log10(final_fit_chi2_p_vals), -np.log10(final_med_p_vals)]
        labels = [r"$\chi^2(1)$", r"$\chi^2(2)$", rf"$\chi^2(k={{{fit_chi2['dist_fit'].args[0]:.2f}}})$", "Gen. gamma", ]
        colors = ["#861A5E", colorlist[0], "r", colorlist[2]]
        axins = axs.inset_axes([.65, .11, .3, .3])
        #axs.text(-.2, .97, r"$\bf{D}$", fontsize=13, transform=axs.transAxes)
        plot_qq(axs, axins, logps, labels, colors=colors, legend_loc = "upper left", thin=True, rasterized=True)
        fig.savefig(Path(f"{output_dir}/neutral_g{num_gens}_d{init_dist}_{output_suffix}{title}_300k_final.pdf"), format="pdf", bbox_inches="tight", dpi=300)
        plt.close(fig)