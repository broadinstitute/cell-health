#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import numpy as np
import pandas as pd

import mygene
from Bio.Seq import Seq


# In[2]:


# Load IDR file
file = "idr0080-screenA-library.csv"

df = pd.read_csv(file)

print(df.shape)
df.head(2)


# In[3]:


# Get ensemble IDs - this will require a manual check
mg = mygene.MyGeneInfo()

result = mg.querymany(
    df.loc[:, "Gene Symbol"].unique().tolist(),
    scopes="symbol,alias",
    species="human",
    fields="entrezgene,symbol,ensembl.gene,",
    as_dataframe=True
)

ensembl_id_df = (
    result
    .sort_values(by="_score", ascending=False)
    .reset_index()
    .drop_duplicates(subset="query")
)


# In[4]:


# Get antisense sequences
# Also, update the comments if the control type is an empty well
antisense_results = []
comments_update = []
for idx, perturbation in df.iterrows():
    sequence = perturbation["Sense Sequence"]
    control_type = perturbation["Control Type"]
    comments = perturbation["Comments"]

    if isinstance(sequence, str):
        my_dna = Seq(sequence)
        antisense = str(my_dna.complement())
    else:
        antisense = np.nan
    
    antisense_results.append(antisense)
    
    if control_type == "empty well":
        comments = f"{comments}; empty well"
    
    comments_update.append(comments)

assert len(antisense_results) == df.shape[0]


# In[5]:


# Drop SMARCB1 (mygene failure)
# Drop query "nan" (mygene false positive)
#{"SMARCB1": "ENSG00000099956"}
ensembl_id_df = ensembl_id_df.dropna(subset=["ensembl.gene"], axis="rows").query("symbol != 'SCN11A'")

ensembl_mapper = dict(zip(ensembl_id_df.loc[:, "query"], ensembl_id_df.loc[:, "ensembl.gene"]))


# In[6]:


# Identify ensembl column
ensemble_column = df.loc[:, "Gene Symbol"].replace(ensembl_mapper)
ensemble_column = [x if str(x).startswith("ENSG") else np.nan for x in ensemble_column]
annotation_build_column = [
    "Ensembl release 101 - August 2020" if str(x).startswith("ENSG") else np.nan for x in ensemble_column
]


# In[7]:


# Update control type
control_type_col = df.loc[:, "Control Type"].replace({"empty well": "no reagent"})


# In[8]:


df.loc[:, "Antisense Sequence"] = antisense_results
df.loc[:, "Gene Identifier"] = ensemble_column
df.loc[:, "Comments"] = comments_update
df.loc[:, "Control Type"] = control_type_col

# Add annotation build
df.loc[:, "Reagent Design Gene Annotation Build"] = annotation_build_column


# In[9]:


# Output updated results
output_file = f"updated_{file}"
df.to_csv(output_file, index=False, sep=",")


# In[10]:


print(df.shape)
df.head(3)

