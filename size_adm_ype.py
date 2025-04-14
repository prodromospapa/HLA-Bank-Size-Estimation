import pandas as pd
# pip install lxml

# administrative units
wiki_page = 'https://en.wikipedia.org/wiki/Regions_of_Greece#List_of_regions'
administrations = pd.read_html(wiki_page)[4]
administrations.set_index(administrations.columns[0],inplace=True)
administrations = administrations[["Region","Population (2024)[2]"]]
administrations.columns = ["Region","Population"]
administrations.drop('(14)',inplace=True) # I don't need it
administrations = administrations.rename_axis("administation")

ype = {1:[1],2:[8,10],3:[3,13],4:[5],5:[11,2],6:[6,7,12,9],7:[4]}

df = pd.DataFrame(columns=["administration_cluster","health_cluster"],index=range(1,14))

df["administration_cluster"] = administrations["Population"].tolist()
for ype,ype_adms in ype.items():
    df["health_cluster"][ype] = df["administration_cluster"][ype_adms].sum()

df.to_pickle("size_adm_ype.pkl")