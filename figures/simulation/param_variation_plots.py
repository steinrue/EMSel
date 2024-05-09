import pickle
import matplotlib.pyplot as plt
from pathlib import Path
from emsel_util import params_dict_to_str, convert_from_abbrevs
from cycler import cycler
from copy import deepcopy
import pandas as pd
import seaborn as sns

###### MODIFY

output_dir = "output"
EM_dir = "EM/param_variation"

###### DO NOT MODIFY

plt.rcParams.update({'font.size': 9,
                     'text.usetex': False,
                     'font.family': 'serif',
                     'font.serif': 'cmr10',
                     'mathtext.fontset': 'cm',
                     'axes.unicode_minus': False,
                     'axes.formatter.use_mathtext': True,})
#graphs = np.where(estimated_s_matrix[:, 1] < 0)[0]
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

coolormap = plt.get_cmap("Dark2")
colorlist = ["#1D6996", *[coolormap(i) for i in [1,0]], colors[3], colors[4]]
plt.rcParams["axes.prop_cycle"] = cycler(color=colorlist)

max_axis = 5

timepoints_list = [2, 5, 11, 35, 101]
samples_list = [6, 20, 50, 100, 200]
knees = [100, 1000, 10000, 100000, 1000000]
hidden_states_list = [100,250,500,1000,2000]
sel_str = .025
onep_sel_types = ["add", "dom", "rec"]
long_sel_names = convert_from_abbrevs(onep_sel_types)
num_gens = 251
init_dist = .25

var_names = ["num_samples", "sample_times", "Ne", "hidden_states"]
label_names = ["Number of samples", "Sampling timepoints", "Population size (Ne)", "Number of HMM hidden states"]
subfig_names = ["A", "B", "C", "D"]
parameter_lists = [samples_list, timepoints_list, knees, hidden_states_list]

knee_sci_ticks = [f"{knees[i]:.0e}" for i in range(len(knees))]

fig, axs = plt.subplots(2, 2, figsize=(6,6), layout="constrained")
for v_i, var_list in enumerate(parameter_lists):
    df_dict = {"s_vals": [], "s_types": [], "var_vals": []}
    v_internal = 0
    for var in var_list:
        fname_dict = {}
        fname_dict["sel_str"] = sel_str
        fname_dict["init_dist"] = init_dist
        fname_dict["num_gens"] = num_gens
        fname_dict[f"{var_names[v_i]}"] = var
        fname = params_dict_to_str(**fname_dict)
        for sel_i, sel_type in enumerate(onep_sel_types):
            hmm_name_dict = deepcopy(fname_dict)
            hmm_name_dict["sel_type"] = sel_type
            hmm_name = params_dict_to_str(**hmm_name_dict)
            hmm_filename = Path(f"{EM_dir}/{hmm_name}_EM.pkl")
            if not hmm_filename.is_file():
                print(f"no file: {hmm_filename}")
                continue
            with open(hmm_filename, "rb") as file:
                hf = pickle.load(file)

            s_vals = hf[f"{sel_type}_run"]["s_final"][1, :]
            df_dict["s_vals"].extend(s_vals.tolist())
            df_dict["s_types"].extend([convert_from_abbrevs(f"{sel_type}", shortall=True)] * s_vals.shape[0])
            df_dict["var_vals"].extend([var] * s_vals.shape[0])
        v_internal += 1
    df = pd.DataFrame(df_dict)
    sns.boxplot(df, x="var_vals", y="s_vals", hue="s_types", dodge=True, ax=axs[v_i//2, v_i%2], native_scale=True,
                log_scale=(True, False), width=.75,
                fliersize=1, boxprops={"lw": .5}, medianprops={"lw": .5}, whiskerprops={"lw": .5}, capprops={"lw": .5},
                flierprops={"alpha": .7}, whis=(2.5, 97.5))
    if not v_i%2:
        axs[v_i//2, v_i%2].set_ylabel(r"$\hat{s}$")
    else:
        axs[v_i//2, v_i%2].set_ylabel("")
    axs[v_i//2, v_i%2].text(-.2, .97, rf"$\bf{subfig_names[v_i]}$", fontsize=13, transform=axs[v_i//2, v_i%2].transAxes)
    axs[v_i//2, v_i%2].axhline(sel_str, color="r", ls="--")
    axs[v_i//2, v_i%2].legend(fontsize=7, labelspacing=.2, handlelength=1.5, handleheight=.5, handletextpad=.4,
                                  borderpad=.2, borderaxespad=.2, markerscale=.25)
    axs[v_i//2, v_i%2].set_ylim([0, 2*sel_str])
    axs[v_i//2, v_i%2].set_xlabel(f"{label_names[v_i]}")
    if v_i != 2:
        axs[v_i//2, v_i%2].set_xticks(parameter_lists[v_i], labels=parameter_lists[v_i])
fig.savefig(f"{output_dir}/param_variation_plot.pdf", format="pdf", bbox_inches="tight")
plt.close(fig)