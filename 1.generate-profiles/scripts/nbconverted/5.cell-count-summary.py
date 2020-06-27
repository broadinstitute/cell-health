#!/usr/bin/env python
# coding: utf-8

# # Describe Cell Count in the Cell Painting Experiment

# In[1]:


import pathlib
import pandas as pd


# In[2]:


cell_count_dir = pathlib.Path("results")

cell_count_df = []
for file in cell_count_dir.iterdir():
    if "cell_count" in str(file):
        plate = file.name.split("_")[0]
        df = pd.read_csv(file, sep="\t").rename({plate: "cell_count"}, axis="columns")
        cell_count_df.append(df)
        
cell_count_df = pd.concat(cell_count_df)

print(cell_count_df.shape)
cell_count_df.head()


# In[3]:


# Total number of cells
cell_count_df.cell_count.sum()


# In[4]:


# Count cells by perturbation
cell_count_df.groupby("pert_name")["cell_count"].sum().sort_values(ascending=False)


# In[5]:


# Cell count by cell line
cell_line_groupby = cell_count_df.groupby("cell_line")
cell_line_groupby["cell_count"].sum()


# In[6]:


cell_line_groupby["cell_count"].hist(bins=30)


# In[7]:


# Also, how many features?
profile_dir = pathlib.Path("data", "profiles")
for file in profile_dir.iterdir():
    plate = file.name
    plate_dir = pathlib.Path(profile_dir, plate)
    if plate != '.DS_Store':
        profile_file = plate_dir / f"{plate}.csv.gz"
        profile_df = pd.read_csv(profile_file)
        print(f"{plate}:")
        print(profile_df.shape)

