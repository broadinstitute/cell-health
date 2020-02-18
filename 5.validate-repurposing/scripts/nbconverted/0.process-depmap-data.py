#!/usr/bin/env python
# coding: utf-8

# # Processing DepMap Data
# 
# The cancer dependency map profiled cell viability in A549 using an orthogonal assay.
# Ultimately, I will determine if the viability predictions are correlated with the observed viability effects.
# 
# Here, I process the DepMap viability data to prepare for analysis input.
# 
# ## Calculating Viability in DepMap
# 
# The text below was taken directly from the Cancer Dependency Map README file `primary-screen-readme.txt`.
# See [`data/`](data/) for more licensing details.
# 
# ```
# The primary PRISM Repurposing dataset contains the results of pooled-cell line chemical-perturbation viability screens for 4,518 compounds screened against 578 or 562 cell lines. It was processed using the following steps:
# 
# - Calculate the median fluorescence intensity (MFI) of each bead-replicate pair.
# - Remove outlier-pools. MFI values are log-transformed and centered to the median logMFI for each cell line on each detection plate. For each well on a plate, the median of the centered values is standardized (using median and MAD) against all other wells in the same position in the same screen. Data from wells with a standardized score greater than 5 or less than -5 is filtered.
# - Remove cell line-plates with low control separation. For each cell line-detection plate combination, the SSMD is calculated between the logMFI of positive and negative control treatments. Cell line-detection plates with an SSMD > -2 are removed.
# - Divide data from each well-detection plate by the median of control barcodes from the same well-detection plate.
# - Calculate log2-fold-change for data relative to negative-control wells from the same cell line on the same detection plate.
# - Positive-control and negative-control treatments are removed.
# - Run ComBat for each treatment condition (compound plate-well combination), considering log-viability values as probes and cell line pool-replicate combinations as batches.
# - Median-collapse each drug-dose-cell line combination
# ```

# In[1]:


import os
import numpy as np
import pandas as pd


# In[2]:


def recode_dose(x, doses, return_level=False):
    closest_index = np.argmin([np.abs(dose - x) for dose in doses])
    if np.isnan(x):
        return 0
    if return_level:
        return closest_index + 1
    else:
        return doses[closest_index]


# In[3]:


data_dir = "data"
cell_line_id = "A549_LUNG"
primary_dose_mapping = [0.04, 0.12, 0.37, 1.11, 3.33, 10, 20]


# ## Load Cell Line Information

# In[4]:


cell_file = os.path.join(data_dir, "primary-screen-cell-line-info.csv")
cell_df = pd.read_csv(cell_file, index_col=0)

print(cell_df.shape)
cell_df.head()


# In[5]:


cell_df = cell_df.query("ccle_name == @cell_line_id")
cell_df


# In[6]:


depmap_id = cell_df.depmap_id.values[0]
depmap_id


# ## Load Treatment Information

# In[7]:


treatment_file = os.path.join(data_dir, "primary-screen-replicate-collapsed-treatment-info.csv")
treatment_df = pd.read_csv(treatment_file)

print(treatment_df.shape)
treatment_df.head()


# ## Load Viability Estimates

# In[8]:


viability_file = os.path.join(data_dir, "primary-screen-replicate-collapsed-logfold-change.csv")
viability_df = pd.read_csv(viability_file, index_col=0)

print(viability_df.shape)
viability_df.head()


# In[9]:


viability_df = (
    pd.DataFrame(
        viability_df
        .loc[depmap_id, :]
    )
    .reset_index()
    .rename(
        {
            "index": "treatment_id",
            depmap_id: cell_line_id
        },
        axis="columns"
    )
)

viability_df.head()


# ## Merge Results

# In[10]:


complete_results_df = (
    viability_df
    .merge(
        treatment_df,
        left_on="treatment_id",
        right_on="column_name",
        how="inner"
    )
)

print(complete_results_df.shape)
complete_results_df.head(3)


# ## Recode Dose Information
# 
# Transform to same scale of Drug Repurposing Hub compound doses.

# In[11]:


complete_results_df = complete_results_df.assign(
    dose_recode=(
        complete_results_df
        .dose
        .apply(
            lambda x: recode_dose(x, primary_dose_mapping, return_level=True)
        )
    )
)

print(complete_results_df.shape)
complete_results_df.head(3)


# ## Output Processed File

# In[12]:


output_file = os.path.join(
    "data",
    "processed",
    "{}_viability_estimates.tsv".format(cell_line_id)
)

complete_results_df.to_csv(output_file, sep='\t', index=False)

