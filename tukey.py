import pandas as pd
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap
from tqdm import tqdm


split_group = "health_cluster"
completeness_list = ["full","8/10","9/10"]

def fix_name(name):
    if type(name) is list:
        if name == ["all"]:
            name = "all"
        elif len(name) == 1:
            name = str(name).replace("]","").replace("[","").replace(" ","")
        else:
            name = str(name).replace(" ","")
    elif type(name) is int:
        name = str(name)
    return name

def tukey(recepient_cluster, compl_file,split_group):
    df = pd.read_csv(f"probabilities/health_cluster/{recepient_cluster}_5{compl_file}.csv",index_col=0)
    if recepient_cluster == "all":
        pop_size = size_adm_ype[split_group].sum()
    else:
        pop_size = size_adm_ype[split_group][eval(recepient_cluster)].sum()
    df = df*pop_size

    df_long = df.reset_index().melt(id_vars="index", var_name="BankSize", value_name="Probability")
    df_long.columns = ["Group", "BankSize", "Probability"]
    # Ensure 'Group' is categorical for the test
    df_long['Group'] = df_long['Group'].astype(str)

    # Apply pairwise Tukey HSD
    tukey_result = pairwise_tukeyhsd(df_long['Probability'], df_long['Group'], alpha=0.05)
    tukey_df = pd.DataFrame(data=tukey_result.summary().data[1:], columns=tukey_result.summary().data[0])
    # Display results
    return tukey_df

size_adm_ype = pd.read_pickle("size_adm_ype.pkl")
if not os.path.exists("tukey"):
    os.makedirs("tukey")

if split_group == "health_cluster":
    groups = [[1,2],[3,4],5,6,7] 

if 'groups' not in globals():
    data = pd.read_pickle("data.pkl")
    groups = np.unique(data[split_group].values).tolist()

groups.append("all")

for completeness, group in tqdm([(c, g) for c in completeness_list for g in groups], desc="Processing completeness levels and groups"):
    if completeness == "9/10":
        compl_file = "_9_10"
    elif completeness == "8/10":
        compl_file = "_8_10"
    else:
        compl_file = ""
    if not os.path.exists(f"tukey/{split_group}_{group}{compl_file}.png"):
        group = fix_name(group)
        result = tukey(group, compl_file, split_group)
        result_pivot = result.pivot(index="group1", columns="group2", values="p-adj")
        sns.heatmap(result_pivot, annot=True, cmap=ListedColormap(['blue', 'red']), cbar=False, linewidths=0.5, vmin=0, vmax=0.05)
        plt.title(f"Tukey HSD Test Results for {split_group} {group} {completeness}")
        plt.tight_layout()
        plt.savefig(f"tukey/{split_group}_{group}{compl_file}.png")
        plt.close()
