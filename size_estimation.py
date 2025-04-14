import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import itertools
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
import re
import pickle
from tqdm import tqdm

import warnings
import matplotlib
warnings.filterwarnings("ignore")
matplotlib.use("Agg")

import gc
import traceback
gc.collect()

import sys
import numexpr as ne

# You need to add arlequin folder in the path 

############################################################################################################
#input
data = pd.read_pickle("data.pkl")

max_threads = 'max'
max_percentage_ram = 0.8
tmp_home = True # if True, the script will create the temporary folders in the home directory
if tmp_home:
    dicts_location = "/local/scratch/s3654095" # location of the dictionaries
split_group = "health_cluster" # group to split the data
completeness_list = ["full","8/10","9/10"] # completeness of the data, full, 8/10, 9/10

max_size = 2*(10**6)
step = max_size//1000
loci_all = [5] # loci to use
if split_group == "health_cluster":
    groups = [[1,2],[3,4],5,6,7] 
#input
############################################################################################################



def xl2arp(data,bank_name,loci):
    #bank_name = bank["bank"].iloc[0]
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

def generate_pairs_with_replacement(iterable):
    pairs = []
    for a, b in itertools.product(iterable, repeat=2):
        pairs.append((a, b))
    return pairs

def hapl_dict(pairs,completeness=None):
    dic = {}
    for _,pair in enumerate(pairs):
        freq = float(pair[0][1]) * float(pair[1][1])
        gen_ = list(zip(pair[0][0], pair[1][0]))
        gen_ = list(map(lambda x: sorted(x), gen_))
        if completeness == "9/10":
            gen_ = [gen_[:4]+[[gen_[4][0]]],gen_[:4]+[[gen_[4][1]]]]
        else:
            gen_ = [gen_]
        for gen_part in gen_:
            gen_part = str(gen_part)
            if gen_part in dic:
                dic[gen_part] += freq
            else:
                dic[gen_part] = freq
    # normalize the frequencies
    if completeness != "9/10":
        total = sum(dic.values())
        for key in dic:
            dic[key] /= total
    return dic

def data2dict(group_combination):
    group = group_combination[0]
    loci = group_combination[1]
    completeness = group_combination[2]
    if completeness == "9/10":
        compl_file = "_9_10"
    elif completeness == "8/10":
        compl_file = "_8_10"
    else:
        compl_file = ""
    if group == "all":
        data_part = data.copy()
    else:
        if type(eval(group)) == list:
            data_part = data[data[split_group].isin(eval(group))]
        elif group.isnumeric():
            data_part = data[data[split_group] == int(group)]

    if loci == 3:
        data_part = data_part.drop(columns = ["C1","C2","DQB1_1","DQB1_2"])
    elif loci == 5:
        data_part = data_part[data_part["loci"]==5]
    if group == "all":
        haplotypes_all = pd.DataFrame(columns=["haplotypes", "frequencies"])
        pending_groups = [g for g in groups if g != ["all"]]
        while pending_groups:
            pending_groups_copy = pending_groups.copy()
            for group_ in pending_groups:
                og_name = group_
                group_ = fix_name(group_)
                if not os.path.exists(f"haplotypes/{split_group}/{group_}_{loci}{compl_file}.pkl") and os.path.exists(f"{dicts_location}/freq_dict/{split_group}/{group_}_{loci}{compl_file}.pkl"):
                    hapl = arp2hapl(data_part,group_,loci,compl_file)
                    hapl.to_pickle(f"haplotypes/{split_group}/{group_}_{loci}{compl_file}.pkl")
                    pending_groups_copy.remove(og_name)
                else:
                    if os.path.exists(f"haplotypes/{split_group}/{group_}_{loci}{compl_file}.pkl"):
                        try:
                            hapl = pd.read_pickle(f"haplotypes/{split_group}/{group_}_{loci}{compl_file}.pkl")
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
    else:
        hapl = arp2hapl(data_part,group,loci,compl_file)
        hapl.to_pickle(f"haplotypes/{split_group}/{group}_{loci}{compl_file}.pkl")

    hapl["haplotypes"] = hapl["haplotypes"].apply(lambda x: x.split("~"))
    if completeness == "8/10":
        hapl["haplotypes"] = hapl["haplotypes"].apply(lambda x: x[:4])
        hapl["haplotypes"] = hapl["haplotypes"].apply(tuple)
        hapl = hapl.groupby('haplotypes', as_index=False).agg({'frequencies': 'sum'})
        hapl["haplotypes"] = hapl["haplotypes"].apply(list)
    hapl_val = hapl.values.tolist()
    pairs = generate_pairs_with_replacement(hapl_val)
    freq = hapl_dict(pairs,completeness)
    with open(f"{dicts_location}/freq_dict/{split_group}/{group}_{loci}{compl_file}.pkl","wb") as f:
        pickle.dump(freq,f)
    if tmp_home:
        os.system(f"rsync -aq {dicts_location}/freq_dict/{split_group}/{group}_{loci}{compl_file}.pkl {curr_dir}/freq_dict/{split_group}/{group}_{loci}{compl_file}.pkl >/dev/null 2>&1")
    return {f"{group}_{loci}{compl_file}":freq}

def arp2hapl(data_part,group,loci,compl_file):
    arp_content = xl2arp(data_part,group,loci)

    with open(f"arlequin_files/{split_group}_{group}_{loci}{compl_file}.arp","w") as f:
        f.write(arp_content)
    while True:
        try:
            os.system(f'bash LaunchArlecore.sh "arlequin_files/{split_group}_{group}_{loci}{compl_file}.arp" > out.log 2> err.log')
            hapl = haplo(f"arlequin_files/{split_group}_{group}_{loci}{compl_file}.res/{split_group}_{group}_{loci}{compl_file}.xml")    
            break
        except Exception:
            continue
    return hapl
    
    

def find_probability(donor_bank,rec_bank,loci,sizes,compl_file,split):
    rec_dict = all_dicts[f"{rec_bank}_{loci}{compl_file}"]

    P_donor = np.array([all_dicts[f"{donor_bank}_{loci}{compl_file}"].get(gt, 0.0) for gt in rec_dict.keys()], dtype=np.float64)[:, np.newaxis]
    P_rec = np.array(list(rec_dict.values()), dtype=np.float64)[:, np.newaxis] 
    del rec_dict 
    if compl_file == "_9_10":
        P_rec = P_rec / P_rec.sum(axis=0)
    result = np.zeros((len(sizes),), dtype=np.float64)
    bionom_donor = 1 - P_donor
    for idx in np.array_split(np.arange(len(sizes)), split):
        sizes_part = np.array(sizes)[idx]
        result[idx] = np.sum(ne.evaluate("P_rec * (bionom_donor ** sizes_part)"), axis=0)
    return 1 - result 

def one_plot_process_banks(rec_name, groups, loci, sizes, completeness):
    global max_threads
    global boots
    if completeness == "9/10":
        compl_file = "_9_10"
    elif completeness == "8/10":
        compl_file = "_8_10"
    else:
        compl_file = ""
    
    # calculate RAM to estimate the size of the data
    rec_dict = all_dicts[f"{rec_name}_{loci}{compl_file}"]
    all_dict_size = deep_getsizeof(all_dicts[f"5_{loci}{compl_file}"])
    P_donor_size = len(rec_dict) * np.dtype(np.float64).itemsize
    P_rec_size = len(rec_dict) * np.dtype(np.float64).itemsize
    P_size = P_donor_size + P_rec_size
    rec_dict_size = deep_getsizeof(rec_dict)
    results_size = len(sizes) * np.dtype(np.float64).itemsize
    ram_usage_no_prob = P_size + rec_dict_size + all_dict_size + results_size

    prob_size = len(rec_dict) * len(sizes) * np.dtype(np.float64).itemsize 
    
    free_ram = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_AVPHYS_PAGES')
    if max_threads == "max" or max_threads >= len(groups):
        free_ram_part = (max_percentage_ram * free_ram) / len(groups)
        max_threads = len(groups)
    else:
        free_ram_part = (max_percentage_ram * free_ram) / max_threads

    split = 1 + (prob_size // (free_ram_part - ram_usage_no_prob))
    
    prob_array = np.zeros((len(groups), len(sizes)), dtype=np.float32)
    
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = {executor.submit(find_probability, fix_name(group), rec_name, loci, sizes, compl_file,split): group for group in groups}
        with tqdm(total=len(groups), desc=f"Creating one-plot for {rec_name} {completeness}") as pbar:
            for future in as_completed(futures):
                try:
                    result = future.result()
                    prob_array[groups.index(futures[future]), :] = result
                    del result
                except Exception:
                    print(f"Error processing {futures[future]}: {traceback.format_exc()}")
                finally:
                        # Explicitly delete the executor to free memory
                        pbar.update(1)
            gc.collect()  # Force garbage collection
    # save it to csv file
    prob_df = pd.DataFrame(prob_array, columns=sizes, index=[fix_name(group) for group in groups])
    prob_df.to_csv(f"{curr_dir}/probabilities/{split_group}/{rec_name}_{loci}{compl_file}.csv")

    # Create a separate figure to avoid conflicts in multithreading
    plt.figure()
    for i, donor_name in enumerate(groups):
        plt.plot(sizes, prob_array[i], label=fix_name(donor_name))
    plt.xlabel("Bank size")
    plt.ylabel("Probability of finding a donor (%)")
    if completeness:
        plt.title(f"Probability of finding at least a donor for {rec_name} at {loci} loci ({completeness})")
    else:
        plt.title(f"Probability of finding at least a donor for {rec_name} at {loci} loci")
    plt.legend()
    plt.grid()
    plt.ylim(0, 1)
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x * 100)}"))
    plt.savefig(f"{curr_dir}/plots/{split_group}/one_plot_{rec_name}_{loci}{compl_file}.png", dpi=300)
    plt.close()


def fix_name(name):
    if type(name) is list:
        if name == ["all"]:
            name = "all"
        elif len(name) == 1:
            name = str(name).replace("]","").replace("[","").replace(" ","")
        else:
            name = str(name).replace(" ","")
    return name

def deep_getsizeof(o, ids=set()):
    if id(o) in ids:
        return 0
    r = sys.getsizeof(o)
    ids.add(id(o))
    if isinstance(o, dict):
        r += sum(deep_getsizeof(k, ids) + deep_getsizeof(v, ids) for k, v in o.items())
    elif isinstance(o, (list, tuple, set, frozenset)):
        r += sum(deep_getsizeof(i, ids) for i in o)
    return r

if not os.path.exists(dicts_location):
    print(f"'{dicts_location}' does not exist, please create it or change the path")
    exit()

curr_dir = os.getcwd()
if split_group not in data.columns:
    print(f'"{split_group}" is not in the data columns choose one of {data.columns}')
    exit()
data = data[data[split_group].notna()]


if 'groups' not in globals():
    groups = np.unique(data[split_group].values).tolist()
groups.append("all")
size_adm_ype = pd.read_pickle("size_adm_ype.pkl")

if tmp_home:
    main_dir = f"{os.path.expanduser('~')}/"
    os.chdir(main_dir)
    if not os.path.exists(f"{curr_dir}/freq_dict"):
        os.makedirs(f"{curr_dir}/freq_dict")
    if not os.path.exists(f"{curr_dir}/freq_dict/{split_group}"):
        os.makedirs(f"{curr_dir}/freq_dict/{split_group}")
else:
    dicts_location = curr_dir

if not os.path.exists(f"{dicts_location}/freq_dict"):
    os.makedirs(f"{dicts_location}/freq_dict")
if not os.path.exists(f"{dicts_location}/freq_dict/{split_group}"):
    os.makedirs(f"{dicts_location}/freq_dict/{split_group}")
if not os.path.exists(f"{curr_dir}/plots"):
    os.makedirs(f"{curr_dir}/plots")
if not os.path.exists(f"{curr_dir}/plots/{split_group}"):
    os.makedirs(f"{curr_dir}/plots/{split_group}")
if not os.path.exists(f"haplotypes/{split_group}"):
    os.makedirs(f"haplotypes/{split_group}")

sizes = [i for i in range(0, max_size, step)] # list of sizes to sample 
groups = [group for group in groups if type(group) == list] + [[group] for group in groups if type(group) != list]
if ["all"] in groups:
    if groups[-1] != ["all"]:
        groups.remove(["all"])
        groups.append(["all"])

dict_combinations = []
for group, locus in itertools.product(groups, loci_all):
    if locus == 5:
        dict_combinations.extend((group, locus, comp) for comp in completeness_list)
    else:
        dict_combinations.append((group, locus, "full"))
all_plots = dict_combinations.copy()
missing_dict_combinations = []
existing_dict_combinations = []
all_dicts = {}
for group, loci, completeness in dict_combinations:
    group = fix_name(group)
    if completeness == "9/10":
        compl_file = "_9_10"
    elif completeness == "8/10":
        compl_file = "_8_10"            
    else:
        compl_file = ""
    if not os.path.exists(f"{dicts_location}/freq_dict/{split_group}/{group}_{loci}{compl_file}.pkl"):        
        if (group, loci, completeness) not in missing_dict_combinations:
            missing_dict_combinations.append((group, loci, completeness))
    else:
        if [group, loci, completeness] not in existing_dict_combinations:
            existing_dict_combinations.append([group, loci, completeness])

# load existin dicts
if existing_dict_combinations:
    for file in tqdm(existing_dict_combinations, desc="Loading dictionaries"):
        group, loci, completeness = file
        if completeness == "9/10":
            compl_file = "_9_10"
        elif completeness == "8/10":
            compl_file = "_8_10"            
        else:
            compl_file = ""
        with open(f"{dicts_location}/freq_dict/{split_group}/{group}_{loci}{compl_file}.pkl","rb") as f:
            try:
                all_dicts[f"{group}_{loci}{compl_file}"] = pickle.load(f)
            except EOFError:
                missing_dict_combinations.append((group, loci, completeness))
                os.remove(f"{dicts_location}/freq_dict/{split_group}/{group}_{loci}{compl_file}.pkl")

# ELB algorithm for missing dicts
if missing_dict_combinations:
    if not os.path.exists(f"arlequin_files"):
        os.makedirs(f"arlequin_files")
    with ThreadPoolExecutor(max_workers=len(missing_dict_combinations)) as executor:
        futures = {executor.submit(data2dict, combination): combination for combination in missing_dict_combinations}
        with tqdm(total=len(missing_dict_combinations), desc="Creating dictionaries") as pbar:
            for future in as_completed(futures):
                try:
                    result = future.result()
                    all_dicts.update(result)
                    del result
                except Exception:
                    print(f"Error processing {futures[future]}: {traceback.format_exc()}")
                finally:
                    # Explicitly delete the executor to free memory
                    pbar.update(1)
        gc.collect()  # Force garbage collection
    os.system(f"rm -r arlequin_files")

if split_group in ["health_cluster","administration_cluster"]:
    adm_ype = True

if not os.path.exists(f"{curr_dir}/probabilities"):
    os.makedirs(f"{curr_dir}/probabilities")
if not os.path.exists(f"{curr_dir}/probabilities/{split_group}"):
    os.makedirs(f"{curr_dir}/probabilities/{split_group}")

if not adm_ype:
    boots = 1

missing_plots = []
for plot in all_plots:
    rec_name, loci, completeness = plot
    rec_name = fix_name(rec_name)
    if completeness == "9/10":
        compl_file = "_9_10"
    elif completeness == "8/10":
        compl_file = "_8_10"            
    else:
        compl_file = ""
    if not os.path.exists(f"{curr_dir}/plots/{split_group}/one_plot_{rec_name}_{loci}{compl_file}.png"):
        one_plot_process_banks(rec_name, groups, loci, sizes, completeness)

os.system(f"rm -r haplotypes")                    
if tmp_home:
    os.chdir(curr_dir)
