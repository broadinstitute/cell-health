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

print(x_consensus_df.shape)
x_consensus_df.head(2)


# In[5]:


file = os.path.join(data_dir, "consensus", "cell_health_modz.tsv.gz")
y_consensus_df = pd.read_csv(file, sep="\t")

print(y_consensus_df.shape)
y_consensus_df.head(2)


# ## Split into Training and Testing

# In[6]:


x_train_df, x_test_df, y_train_df, y_test_df = train_test_split(
    x_consensus_df,
    y_consensus_df,
    test_size=test_proportion,
    stratify=y_consensus_df.Metadata_cell_line,
    random_state=42
)


# In[7]:


print(x_train_df.shape)
print(x_test_df.shape)


# In[8]:


file = os.path.join("data", "x_train.tsv.gz")
x_train_df.to_csv(file, sep="\t", index=False)

file = os.path.join("data", "y_train.tsv.gz")
y_train_df.to_csv(file, sep="\t", index=False)

file = os.path.join("data", "x_test.tsv.gz")
x_test_df.to_csv(file, sep="\t", index=False)

file = os.path.join("data", "y_test.tsv.gz")
y_test_df.to_csv(file, sep="\t", index=False)

