#!/usr/bin/env python
# coding: utf-8

# # Build Table Describing Perturbations Used in this Experiment

# In[1]:


import os
import pandas as pd


# In[2]:


cell_health_dir = os.path.join("data", "cell_health_guides")
cell_painting_dir = os.path.join("data", "metadata", "platemap")


# ## Load and Process Cell Health Data

# In[3]:


cell_health_cols = [
    "gene_name",
    "pert_name",
    "Clone ID",
    "Pathway/Cellular Process",
]

cell_painting_cols = [
    "gene_name",
    "pert_name",
    "broad_sample",
]

cell_lines = ["A549", "ES2", "HCC44"]


# In[4]:


cell_health_file = os.path.join(cell_health_dir, "HCI_BenchmarksV1-2_platemap.csv")

cell_health_platemap_df = pd.read_csv(cell_health_file, sep=",")
cell_health_platemap_df.loc[:, "gene_name"] = (
    cell_health_platemap_df.loc[:, "gene_name"].str.replace("_", "-")
)
cell_health_platemap_df.loc[:, "pert_name"] = (
    cell_health_platemap_df.loc[:, "pert_name"].str.replace("_", "-")
)
print(cell_health_platemap_df.shape)
cell_health_platemap_df.head(2)


# In[5]:


cell_health_annotation_file = os.path.join(cell_health_dir, "HCI_benchmarks_sgRNA_annotation.xlsx")

cell_health_annotation_df = pd.read_excel(cell_health_annotation_file, sep=",")
cell_health_annotation_df.loc[:, "guide id"] = (
    cell_health_annotation_df.loc[:, "guide id"].str.replace("_", "-")
)
cell_health_annotation_df.loc[:, "gene symbol"] = (
    cell_health_annotation_df.loc[:, "gene symbol"].str.replace("_", "-")
)

cell_health_annotation_df.loc[
    cell_health_annotation_df.loc[:, "gene symbol"].isna(),
    "gene symbol"
] = (
    cell_health_annotation_df.loc[:, "guide id"]
)

print(cell_health_annotation_df.shape)
cell_health_annotation_df.head(2)


# In[6]:


cell_health_df = (
    cell_health_platemap_df.merge(
        cell_health_annotation_df,
        left_on=["gene_name", "pert_name"],
        right_on=["gene symbol", "guide id"],
        how="outer"
    )
    .loc[:, cell_health_cols]
    .drop_duplicates()
    .sort_values(by="gene_name")
    .reset_index(drop=True)
    .assign(cell_health_data=True)
    .dropna(subset=["gene_name"])
)

print(cell_health_df.shape)
cell_health_df.head(2)


# ## Load and Process Cell Painting Data

# In[7]:


cp_platemap_file = os.path.join(cell_painting_dir, "DEPENDENCIES1_A549.csv")
cp_platemap_df = (
    pd.read_csv(cp_platemap_file, sep=',')
    .loc[:, cell_painting_cols]
    .drop_duplicates()
    .assign(cell_painting_data=True)
)

print(cp_platemap_df.shape)
cp_platemap_df.head()


# In[8]:


full_df = (
    cell_health_df.merge(
        cp_platemap_df,
        left_on=["gene_name", "pert_name", "Clone ID"],
        right_on=["gene_name", "pert_name", "broad_sample"],
        how="outer"
    )
    .rename(
        {
            "Pathway/Cellular Process": "cell_process"
        },
        axis="columns"
    )
    .loc[:,
         [
             "gene_name",
             "pert_name",
             "broad_sample",
             "cell_process",
             "cell_health_data",
             "cell_painting_data"
         ]
        ]
    .sort_values(by=["gene_name", "pert_name"])
)

full_df.loc[full_df.cell_health_data.isna(), "cell_health_data"] = False
full_df.loc[full_df.cell_painting_data.isna(), "cell_painting_data"] = False

print(full_df.shape)
full_df.head()


# ## Output Supplementary Table

# In[9]:


output_file = os.path.join("tables", "supplementary_table_1_perturbation_details.tsv")
full_df.to_csv(output_file, sep='\t', index=False)

