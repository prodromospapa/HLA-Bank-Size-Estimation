import pandas as pd
import numpy as np
import re
import xml.etree.ElementTree as ET
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
import itertools
import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import traceback

################################
data = pd.read_pickle("data.pkl")
max_threads = 'max'
threshold = 0 # threshold for haplotype frequencies
tmp_home = True # if True, the script will create the temporary folders in the home directory
split_group = "health_cluster"
completeness_list = ["full","_8_10"] # completeness of the data, full, 8/10, 9/10
loci_all = [5]
if split_group == "health_cluster":
    groups = [[1,2],[3,4],5,6,7] 
################################

def xl2arp(data,bank_name,loci):
    arp_content = f'''[Profile]

    Title="Genetic Data"
    NbSamples=1
    GenotypicData=1
    GameticPhase=0
    DataType=STANDARD 
    LocusSeparator=TAB 
    RecessiveData=0
    MissingData='?'

[Data]

[[Samples]]
'''
    counter = 0
    if loci == 3:
        loci_l = [['A1','B1','DRB1_1'],['A2','B2','DRB1_2']]
    else:
        loci_l = [['A1','B1','C1','DRB1_1','DQB1_1'],['A2','B2','C2','DRB1_2','DQB1_2']]
    arp_content += f'''SampleName="{bank_name}"
    SampleSize={len(data)}
    SampleData= {{
'''
    for _, row in data.iterrows():
        counter += 1
        id = row['ID']
        arp_content += f"{id}\t" + f"1\t" + '\t'.join(row[loci_l[0]]) + "\n" + "\t" + '\t'.join(row[loci_l[1]]) + "\n"
    arp_content += "}\n\n"
    arp_content += f'''[[Structure]]
    StructureName="Simulated data"
    NbGroups=1
    Group={{
    "{bank_name}"
    }}'''
    return arp_content

def haplo(file):
    tree = ET.parse(file)
    root = tree.getroot()
    hapl = root[20].text.split("\n\n")[5].split("\n")[:-1]
    lista = [["~".join([re.match(r"^.*\*\d+:\d+", allele).group() if re.match(r"^.*\*\d+:\d+", allele) else allele for allele in sample.split()[4:]])]+sample.split()[2:3] for sample in hapl]
    df = pd.DataFrame(lista,columns=["haplotypes","frequencies"])
    df["frequencies"] = df["frequencies"].astype(np.float32)
    df = df.groupby("haplotypes")["frequencies"].sum().reset_index()
    df["frequencies"] = df["frequencies"] / df["frequencies"].sum()
    return df

def arp2hapl(data_part,group,loci,compl_file):
    if loci == 3:
        data_part = data_part.drop(columns = ["C1","C2","DQB1_1","DQB1_2"])
    elif loci == 5:
        data_part = data_part[data_part["loci"]==5]

    group = fix_name(group)
    if group != "all":
        if type(eval(group)) == list:
            data_part = data_part[data_part[split_group].isin(eval(group))]
        elif group.isnumeric():
            data_part = data_part[data_part[split_group] == int(group)]
        if not os.path.exists(f"haplotypes/{split_group}_{group}_{loci}{compl_file}.pkl"):
            arp_content = xl2arp(data_part,group,loci)
            with open(f"arlequin_files/{split_group}_{group}_{loci}{compl_file}.arp","w") as f:
                f.write(arp_content)
            while True:
                try:
                    os.system(f'bash LaunchArlecore.sh "arlequin_files/{split_group}_{group}_{loci}{compl_file}.arp" > out.log 2> err.log')
                    hapl = haplo(f"arlequin_files/{split_group}_{group}_{loci}{compl_file}.res/{split_group}_{group}_{loci}{compl_file}.xml")  
                    hapl.to_pickle(f"haplotypes/{split_group}_{group}_{loci}{compl_file}.pkl")  
                    break
                except Exception:
                    #print(f"Error processing {group}, {loci}, {compl_file}: {traceback.format_exc()}")
                    continue
        else:
            hapl = pd.read_pickle(f"haplotypes/{split_group}_{group}_{loci}{compl_file}.pkl")
    else:
        pd.set_option('future.no_silent_downcasting', True)
        haplotypes_all = pd.DataFrame(columns=["haplotypes", "frequencies"])
        pending_groups = [g for g in groups if g != ["all"]]
        while pending_groups:
            pending_groups_copy = pending_groups.copy()
            for group_ in pending_groups:
                og_name = group_
                group_ = fix_name(group_)
                if os.path.exists(f"haplotypes/{split_group}_{group_}_{loci}{compl_file}.pkl"):
                    try:
                        hapl = pd.read_pickle(f"haplotypes/{split_group}_{group_}_{loci}{compl_file}.pkl")
                        pending_groups_copy.remove(og_name)
                    except EOFError:
                        continue
                else:
                    continue    
                pop_size = size_adm_ype[split_group][eval(group_)].sum()
                hapl["frequencies"] = hapl["frequencies"] * (pop_size / size_adm_ype[split_group].sum())
                result = pd.merge(haplotypes_all, hapl, on='haplotypes', how='outer', suffixes=('_all', '_haplotypes'))
                result['frequencies'] = result['frequencies_all'].fillna(0) + result['frequencies_haplotypes'].fillna(0)
                haplotypes_all = result[['haplotypes', 'frequencies']]
            pending_groups = pending_groups_copy.copy()
        haplotypes_all["frequencies"] = haplotypes_all["frequencies"] / haplotypes_all["frequencies"].sum()
        hapl = haplotypes_all
    if compl_file == "_8_10":
        hapl["haplotypes"] = hapl["haplotypes"].apply(lambda x: x.split("~"))
        hapl["haplotypes"] = hapl["haplotypes"].apply(lambda x: x[:4])
        hapl["haplotypes"] = hapl["haplotypes"].apply(tuple)
        hapl = hapl.groupby('haplotypes', as_index=False).agg({'frequencies': 'sum'})
        hapl["haplotypes"] = hapl["haplotypes"].apply(list)
        hapl["haplotypes"] = hapl["haplotypes"].apply(lambda x: "~".join(x))
    if threshold > 0:
        hapl = hapl[hapl["frequencies"] >= threshold] 
        hapl["frequencies"] = hapl["frequencies"] / hapl["frequencies"].sum()
    return hapl

def fix_name(name):
    if type(name) is list:
        if name == ["all"]:
            name = "all"
        elif len(name) == 1:
            name = str(name).replace("]","").replace("[","").replace(" ","")
        else:
            name = str(name).replace(" ","")
    return name


# check if threshold necessary
def sorensen(hapl_rec,hapl_donor):
    merged_df = pd.merge(hapl_rec, hapl_donor, on='haplotypes', how='inner')
    inter_sum = (merged_df['frequencies_x'].astype(float) - merged_df['frequencies_y'].astype(float)).abs().sum()
    result = (2*(len(merged_df)-inter_sum)) / (len(hapl_rec) + len(hapl_donor))
    return result

size_adm_ype = pd.read_pickle("size_adm_ype.pkl")

curr_dir = os.getcwd()
if tmp_home:
    main_dir = f"{os.path.expanduser('~')}/"
    os.chdir(main_dir)

data = data[data[split_group].notna()]
if 'groups' not in globals():
    groups = np.unique(data[split_group].values).tolist()
groups.append("all")
groups = [group for group in groups if type(group) == list] + [[group] for group in groups if type(group) != list]
if ["all"] in groups:
    if groups[-1] != ["all"]:
        groups.remove(["all"])
        groups.append(["all"])


if not os.path.exists("arlequin_files"):
    os.mkdir("arlequin_files")
if not os.path.exists(f"{curr_dir}/sorensen"):
    os.mkdir(f"{curr_dir}/sorensen")
if not os.path.exists(f"{curr_dir}/sorensen/{split_group}"):
    os.mkdir(f"{curr_dir}/sorensen/{split_group}")
if not os.path.exists("haplotypes"):
    os.mkdir("haplotypes")

haplotypes = []
for group, locus in itertools.product(groups, loci_all):
    if locus == 5:
        haplotypes.extend((group, locus, comp) for comp in completeness_list)
    else:
        haplotypes.append((group, locus, "full"))

haplo_dict = {locus: {compl: {fix_name(group): None for group in groups} for compl in completeness_list} for locus in loci_all}

if max_threads == "max":
        max_threads = len(haplotypes)
max_threads = min(max_threads, os.cpu_count()//2)
with ThreadPoolExecutor(max_workers=max_threads) as executor:
    futures = {executor.submit(arp2hapl, data, group, locus, compl): (group, locus, compl) for group, locus, compl in haplotypes}
    with tqdm.tqdm(total=len(futures)) as pbar:
        for future in as_completed(futures):
            group, locus, compl = futures[future]
            try:
                hapl = future.result()
                if compl in completeness_list:
                    haplo_dict[locus][compl][fix_name(group)] = hapl
            except Exception as e:
                print(f"Error processing {group}, {locus}, {compl}: {traceback.format_exc()}")
            pbar.update(1)

for loci in loci_all:
    for compl in completeness_list:
        table = pd.DataFrame(columns=haplo_dict[loci][compl].keys(),
                            index=list(haplo_dict[loci][compl].keys())[::-1],dtype=float)
        hapl_combinations = list(itertools.combinations(haplo_dict[loci][compl].keys(), 2))
        for rec, donor in hapl_combinations:
            table.loc[rec,donor] = sorensen(haplo_dict[loci][compl][donor], haplo_dict[loci][compl][rec])
        table = table * 100
        table.to_csv(f"{curr_dir}/sorensen/{split_group}/sorensen_{loci}_{compl}.csv")
        sns.heatmap(table, annot=True, fmt=".2f", cbar=False, mask=table.isnull(), linewidths=0.5,vmin=0, vmax=100)
        if compl == "_8_10":
            compl_title = "8/10"
        else:
            compl_title = "full"
        plt.title(f"Sorensen % similarity {split_group} for {loci} loci ({compl_title})")
        plt.savefig(f"{curr_dir}/sorensen/{split_group}/sorensen_{loci}_{compl}.png")
        plt.close()
os.system("rm -r arlequin_files")
#os.system("rm -r haplotypes")
if tmp_home:
    os.chdir(curr_dir)
