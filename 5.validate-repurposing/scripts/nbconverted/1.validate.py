#!/usr/bin/env python
# coding: utf-8

# # Observe Viability Predictions against DepMap Estimates

# In[1]:


import os
import pathlib
import pandas as pd
from scipy import stats

import plotnine as gg


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


cell_line_id = "A549_LUNG"
cell_health_model = "cell_health_modz_target_cc_cc_n_objects"


# ## Load Data
# 
# ### Cancer Dependency Map Data

# In[4]:


depmap_file = os.path.join(
    "data",
    "processed",
    "{}_viability_estimates.tsv".format(cell_line_id)
)

depmap_df = (
    pd.read_csv(depmap_file, sep='\t')
    .rename(
        {
            cell_line_id: "depmap_viability"
        },
        axis="columns"
    )
)

print(depmap_df.shape)
depmap_df.head(3)


# ### Cell Health Predictions on Drug Repurposing Hub

# In[5]:


focus_columns = [
    "Metadata_Plate_Map_Name",
    "Metadata_broad_sample",
    "Metadata_dose_recode",
    "Metadata_mmoles_per_liter",
    "Metadata_pert_well",
    cell_health_model
]


# In[6]:


cell_health_file = os.path.join(
    "..",
    "4.apply",
    "data",
    "repurposing_transformed_real_models_modz.tsv.gz"
)

cell_health_df = (
    pd.read_csv(cell_health_file, sep='\t')
    .rename(
        {
            cell_health_model: "cell_health_viability"
        },
        axis="columns"
    )
)

print(cell_health_df.shape)
cell_health_df.head(3)


# ## Observe how different dose recodings are

# In[7]:


depmap_df.plot(x="dose", y="dose_recode", kind="scatter")


# In[8]:


cell_health_df.plot(x="Metadata_mmoles_per_liter", y="Metadata_dose_recode", kind="scatter")


# ## Merge Together

# In[9]:


full_df = (
    cell_health_df
    .merge(
        depmap_df,
        left_on=["Metadata_broad_sample", "Metadata_dose_recode"],
        right_on=["broad_id", "dose_recode"],
        how="inner"
    )
    .dropna(subset=["depmap_viability"])
)

output_file = pathlib.Path("results", "depmap_viability_validation.tsv.gz")
full_df.to_csv(output_file, sep="\t", index=False)

print(full_df.shape)
full_df.head(3)


# In[10]:


len(full_df.Metadata_broad_sample.unique())


# In[11]:


dose_differences_gg = (
    gg.ggplot(full_df, gg.aes(x="dose", y="Metadata_mmoles_per_liter")) +
    gg.geom_point(size = 0.5, alpha = 0.6) +
    gg.theme_bw() +
    gg.geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    gg.xlab("PRISM Dose") +
    gg.ylab("Dependency Map Dose") +
    gg.coord_fixed()
)

output_file = os.path.join("figures", "dose_differences.png")
dose_differences_gg.save(output_file, dpi = 400, height = 3, width = 3)

dose_differences_gg


# ## Obtain Results

# In[12]:


spearman_cor = stats.spearmanr(full_df.cell_health_viability, full_df.depmap_viability)
spearman_cor = pd.DataFrame(spearman_cor, index=["stat", "p"]).transpose()
spearman_cor


# In[13]:


result_text = "Spearman = {0:.2f}\np = {1:.2E}".format(
    spearman_cor.stat[0],
    spearman_cor.p[0]
)

result_text


# In[14]:


viability_gg = (
    gg.ggplot(full_df, gg.aes(x="cell_health_viability", y="depmap_viability")) +
    gg.geom_point(size = 0.5, alpha = 0.3) +
    gg.theme_bw() +
    gg.geom_smooth() +
    gg.annotate("text", label = result_text, x = -2, y = -8.5) +
    gg.xlab("Cell Health Model Predictions (Number of Objects)") +
    gg.ylab("Dependency Map Viability") +
    gg.coord_fixed()
)

output_file = os.path.join("figures", "viability_results.png")
viability_gg.save(output_file, dpi = 400, height = 3, width = 5)

viability_gg

