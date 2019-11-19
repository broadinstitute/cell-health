#!/usr/bin/env python
# coding: utf-8

# # Stratify Data into Training and Testing Sets
# 
# **Gregory Way, 2019**
# 
# Split the input data into training and testing sets balanced by cell line.

# In[1]:


import os
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from pycytominer.get_na_columns import get_na_columns


# In[2]:


np.random.seed(123)


# In[3]:


test_proportion = 0.15
data_dir = os.path.join("..", "1.generate-profiles", "data")


# ## Load Data

# In[4]:


file = os.path.join(data_dir, "consensus", "cell_painting_modz.tsv.gz")
x_consensus_df = pd.read_csv(file, sep="\t")

num_original_features = x_consensus_df.shape[1]

print(x_consensus_df.shape)
x_consensus_df.head(2)


# In[5]:


file = os.path.join(data_dir, "consensus", "cell_health_modz.tsv.gz")
y_consensus_df = pd.read_csv(file, sep="\t")

print(y_consensus_df.shape)
y_consensus_df.head(2)


# ## Subset Features into those acquired in Repurposing Data

# In[6]:


# Note, these files are not yet public!
repurposing_project_id = "2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad"
example_plate = "SQ00015058"

repurposing_profile_dir = os.path.join(
    "/Users",
    "gway",
    "work",
    "projects",
    repurposing_project_id,
    "workspace",
    "software",
    repurposing_project_id,
    "subsampling",
    "data",
    "profiles"
)

plate_dir = os.path.join(repurposing_profile_dir, example_plate, "n_all")
example_plate_file = os.path.join(plate_dir, "{}_subsample_all_normalized.csv".format(example_plate))
repurposing_df = pd.read_csv(example_plate_file)

print(repurposing_df.shape)
repurposing_df.head()


# In[7]:


cp_features = set(repurposing_df.columns[~repurposing_df.columns.str.startswith("Metadata")])
cp_features = sorted(
    list(
        cp_features
        .intersection(
            set(
                x_consensus_df.columns[~x_consensus_df.columns.str.startswith("Metadata")]
            )
        )
    )
)

len(cp_features)


# In[8]:


meta_cols = x_consensus_df.columns[x_consensus_df.columns.str.startswith("Metadata")].tolist()
x_consensus_df = x_consensus_df.loc[:, meta_cols + cp_features]
num_subset_features = x_consensus_df.shape[1]

print(x_consensus_df.shape)
x_consensus_df.head()


# In[9]:


print("subsetting by repurposing features removed {} features".format(num_original_features - num_subset_features))


# ## Split into Training and Testing

# In[10]:


x_train_df, x_test_df, y_train_df, y_test_df = train_test_split(
    x_consensus_df,
    y_consensus_df,
    test_size=test_proportion,
    stratify=y_consensus_df.Metadata_cell_line,
    random_state=42
)


# In[11]:


print(x_train_df.shape)
print(x_test_df.shape)


# In[12]:


file = os.path.join("data", "x_train.tsv.gz")
x_train_df.to_csv(file, sep="\t", index=False)

file = os.path.join("data", "y_train.tsv.gz")
y_train_df.to_csv(file, sep="\t", index=False)

file = os.path.join("data", "x_test.tsv.gz")
x_test_df.to_csv(file, sep="\t", index=False)

file = os.path.join("data", "y_test.tsv.gz")
y_test_df.to_csv(file, sep="\t", index=False)

