#!/usr/bin/env python
# coding: utf-8

# # Analyzing Audits
# 
# **Predicting Cell Health**
# 
# **Gregory Way, 2019**
# 
# Detecting replicate correlation across treatments, wells, and cell lines.

# In[1]:


import os
import pandas as pd

import matplotlib.pyplot as plt
import plotnine as gg


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


# Load Audits
# (both detailed treatments and null threshold summaries)
data_dir = "data"

audit_df_list = []
audit_summary_list = []
for file in os.listdir(data_dir):
    full_file = os.path.join(data_dir, file)
    df = pd.read_csv(full_file)
    if 'detailed' in file:
        audit_df_list.append(df)
    if 'audit.csv' in file:
        audit_summary_list.append(df)
        
audit_df = pd.concat(audit_df_list)
summary_audit_df = pd.concat(audit_summary_list).reset_index(drop=True)

print(audit_df.shape)
audit_df.head()


# ## Interpretation of `correlation`
# 
# This represents the median pairwise Pearson correlation of replicates across genes, perturbations, and cell lines.
# The data here represent correlations across cell painting profiles.

# In[4]:


# Process the summary audits
summary_audit_df = (
    pd.concat(
        [summary_audit_df,
         summary_audit_df.plate_map_name.str.split('_', expand=True).iloc[:, 1]],
        axis='columns'
    )
    .rename(columns={1: 'Metadata_cell_line'})
)

summary_audit_df


# ## Interpretation of `null_threshold`
# 
# The `null_threshold` column above is generated in [broadinstitute/cytominer_scripts/audit.R](https://github.com/broadinstitute/cytominer_scripts/blob/master/audit.R).
# 
# The metric represents a random sampling of different samples and combined, over many iterations, to represent a null distribution of correlation.

# In[5]:


# Summarize the audits using median correlation per gene and cell line
median_collapsed_correlation_df = (
    audit_df
    .groupby(['Metadata_gene_name', 'Metadata_pert_name', 'Metadata_cell_line'])
    ['correlation']
    .median()
    .reset_index()
    .groupby(['Metadata_gene_name', 'Metadata_cell_line'])
    ['correlation']
    .median()
    .reset_index()
    .rename(columns={'correlation': 'median_correlation'})
)

median_collapsed_correlation_df.head()


# ## Interpretation of `median_correlation`
# 
# This value represents the median correlation across all **guides** targetting the given gene for the given cell line.
# 
# Therefore, the `median_correlation` is a median (guides) of medians (replicates of single guides) of medians (profiles of individual cells in a single well).

# In[6]:


# Merge summarized data into audit dataframe
merge_cols = ['Metadata_gene_name', 'Metadata_cell_line']

audit_df = (
    audit_df
    .merge(
        median_collapsed_correlation_df,
        left_on=merge_cols,
        right_on=merge_cols
    )
)

# Determine how many CRISPRs are targetting the gene
audit_df = (
    audit_df
    .assign(
        guide_index=(
            audit_df
            .groupby(merge_cols)
            .cumcount()+1
        )
        .tolist()
    )
    .query("Metadata_gene_name != 'Chr2'")
    .sort_values(by='median_correlation',
                 ascending=False)
)

audit_df.guide_index = audit_df.guide_index.astype('str')
audit_df.head()


# In[7]:


# Get audit_df ready to plot
gene_name_cat = pd.CategoricalDtype(categories=audit_df.Metadata_gene_name.unique().tolist())
audit_df = audit_df.assign(
    gene_plot_id=audit_df.Metadata_gene_name.astype(str).astype(gene_name_cat)
)

# Merge with the summary to add in cell line specific lines
audit_df = (
    audit_df
    .merge(summary_audit_df,
           left_on='Metadata_cell_line',
           right_on='Metadata_cell_line')
)

audit_df.head()


# In[8]:


audit_df.head()


# ## Plot Figures

# In[9]:


gg.options.figure_size=(6.4, 4.8)

# Make sure to drop duplicates of redundant gene, perturbation, and cell line columns
# Not removing replicates will put more weight on genes with more measurements

cor_density_gg = (
    gg.ggplot(audit_df.drop_duplicates(["Metadata_gene_name", "Metadata_pert_name", "Metadata_cell_line"]),
              gg.aes(x="median_correlation")) + \
        gg.geom_density(gg.aes(fill="Metadata_cell_line"),
                        alpha=0.4) + \
        gg.theme_bw() + \
        gg.xlim([-0.5, 1]) + \
        gg.xlab("Median Guide Correlation") + \
        gg.ylab("Density") + \
        gg.scale_fill_manual(name="Cell Line",
                             values=["#1b9e77", "#d95f02", "#7570b3"])
)

cor_density_gg


# In[10]:


file = os.path.join("figures", "median-guide-correlation-density")
for extension in ['.png', '.pdf']:
    gg.ggsave(cor_density_gg,
              filename='{}{}'.format(file, extension),
              dpi=300,
              height=4.8,
              width=6.4,
              units='in')


# In[11]:


gg.options.figure_size=(8, 15)

gene_correlation_gg = (
    gg.ggplot(audit_df,
              gg.aes(x="correlation",
                     y="gene_plot_id")) + \
        gg.geom_jitter(size=0.3) + \
        gg.geom_boxplot() + \
        gg.xlab("Pearson Correlation") + \
        gg.ylab("Target") + \
        gg.geom_vline(gg.aes(xintercept="null_threshold"),
                      linetype='dashed',
                      color='red') + \
        gg.xlim([-0.5, 1]) + \
        gg.facet_grid("~Metadata_cell_line") + \
        gg.theme_bw()
)
    
gene_correlation_gg


# In[12]:


file = os.path.join("figures", "target-guide-correlation")
for extension in ['.png', '.pdf']:
    gg.ggsave(gene_correlation_gg,
              filename='{}{}'.format(file, extension),
              dpi=300,
              height=15,
              width=8,
              units='in')

