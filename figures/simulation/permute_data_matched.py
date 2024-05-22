import numpy as np
from pathlib import Path
from tqdm import tqdm


####### MODIFY

data_dir = "data"

###### DO NOT MODIFY
df_path = Path(f"{data_dir}/neutral_g125_dal_special_data.csv")
df = np.loadtxt(df_path, delimiter="\t")
final_data = df[:, 2::3].astype(int)
num_samples = df[:, 1::3].astype(int)
sample_times = (df[:, ::3]).astype(int)
for i in tqdm(range(100)):
    rng = np.random.default_rng(i)
    new_final_data = np.zeros_like(final_data)
    alt_df_path = Path(f"{data_dir}/neutral_g125_dal_special_perm{i}_data.csv")
    for s_i in np.arange(sample_times.shape[0]):
        samples = np.zeros(np.sum(num_samples[s_i, :]))
        new_sample_times = np.zeros_like(samples)
        samples_counter = 0
        for st_i in np.arange(num_samples.shape[1]):
            samples_time_i = num_samples[s_i, st_i]
            hits_time_i = final_data[s_i, st_i]
            samples[samples_counter:samples_counter+hits_time_i] += 1
            new_sample_times[samples_counter:samples_counter+samples_time_i] = sample_times[s_i, st_i]
            samples_counter += samples_time_i
        rng.shuffle(new_sample_times)
        for st_i in np.arange(num_samples.shape[1]):
            new_final_data[s_i, st_i] = np.sum(samples[new_sample_times == sample_times[s_i, st_i]])
    df[:, 2::3] = new_final_data.astype(int)
    np.savetxt(alt_df_path, df, delimiter="\t", fmt="%d",
        header="Each row = one replicate; each set of three columns = (sampling time, total samples, derived alleles)")