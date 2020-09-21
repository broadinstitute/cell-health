#!/usr/bin/env python
# coding: utf-8

# ## Apply UMAP to visualize relationships between the CRISPR perturbations

# In[9]:


import umap
import pathlib
import numpy as np
import pandas as pd
import plotnine as gg

from pycytominer.cyto_utils import infer_cp_features


# In[2]:


np.random.seed(123)


# In[3]:


# Set constants and file names
consensus = "modz"

data_dir = pathlib.Path("..", "1.generate-profiles", "data", "consensus")
results_dir = pathlib.Path("results")

profile_file = pathlib.Path(data_dir, f"cell_painting_{consensus}.tsv.gz")
pred_file = pathlib.Path(results_dir, f"all_model_predictions_{consensus}.tsv")
output_file = pathlib.Path(results_dir, f"profile_umap_with_predictions_{consensus}.tsv")


# In[4]:


# Load profile data
df = (
    pd.read_csv(profile_file, sep="\t")
    .sort_values(by="Metadata_profile_id")
    .reset_index(drop=True)
)

cp_features = infer_cp_features(df)

print(df.shape)
df.head()


# In[5]:


# Load cell health model predictions
pred_df = (
    pd.read_csv(pred_file, sep="\t")
    .sort_values(by="Metadata_profile_id")
    .reset_index(drop=True)
)

print(pred_df.shape)
pred_df.head()


# In[6]:


# Ensure data and predictions are aligned
assert df.Metadata_profile_id.tolist() == pred_df.Metadata_profile_id.tolist()


# In[7]:


# Apply UMAP
reducer = umap.UMAP(random_state=1234, n_components=2)

predict_embedding_df = pd.DataFrame(
    reducer.fit_transform(df.loc[:, cp_features]),
    columns=["umap_x", "umap_y"]
)

predict_embedding_df = (
    predict_embedding_df
    .merge(
        pred_df,
        left_index=True,
        right_index=True
    )
)

print(predict_embedding_df.shape)
predict_embedding_df.head()


# In[8]:


# Output to file
predict_embedding_df.to_csv(output_file, sep="\t", index=False)

