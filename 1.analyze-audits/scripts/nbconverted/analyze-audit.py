#!/usr/bin/env python
# coding: utf-8

# # Analyzing Audits
# 
# **Predicting Cell Health**
# 
# **Gregory Way, 2019**
# 
# Detecting replicate correlation across replicates and CRISPR guides targeting the same gene.
# 
# Audits were generated using `generate_audits.sh`, which calls [`broadinstitute/cytominer_scripts/audit.R`](https://github.com/broadinstitute/cytominer_scripts/blob/master/audit.R).

# In[1]:


import os
import pandas as pd

import matplotlib.pyplot as plt
import plotnine as gg

from scripts.audit_utils import get_confidence


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[3]:


# Load Audits
# (both detailed treatments and null threshold summaries)
data_dir = os.path.join("..", "..", "..", "scratch", "CRISPR_PILOT_B1", "audit")

audit_df_list = []
audit_summary_list = []
audit_guide_df_list = []

# Load all audit data and allocate to proper data frame
for file in os.listdir(data_dir):
    full_file = os.path.join(data_dir, file)
    df = pd.read_csv(full_file)

    if 'detailed.csv' in file:
        audit_df_list.append(df)
    if 'audit.csv' in file:
        audit_summary_list.append(df)
    if 'detailed_guide.csv' in file:
        audit_guide_df_list.append(df)

# Concatenate all detailed audits
audit_df = (
    pd.concat(audit_df_list)
    .sort_values(by=['Metadata_gene_name',
                     'Metadata_pert_name',
                     'Metadata_cell_line'])
    .reset_index(drop=True)
)

# Concatenate the summary audits
summary_audit_df = (
    pd.concat(audit_summary_list)
    .reset_index(drop=True)
)

# Concatenate replicate guide audit
audit_guide_df = (
    pd.concat(audit_guide_df_list)
    .sort_values(by=['Metadata_cell_line',
                     'Metadata_gene_name'])
    .reset_index(drop=True)
)


# In[4]:


# How many genes have CRISPR guides?
print(len(audit_df.Metadata_gene_name.unique()))

# How many guides per gene and per cell line
pd.crosstab(audit_df.Metadata_cell_line,
            audit_df.Metadata_gene_name).T.head(5)


# ## Interpretation of `correlation`
# 
# The audit process (generated using [broadinstitute/cytominer_scripts/audit.R](https://github.com/broadinstitute/cytominer_scripts/blob/master/audit.R)) identifies median pairwise correlations between all replicates within plate maps, genes, and guides.
# 
# Therefore, the `correlation` column represents the pairwise correlation across replicates **within batch**.

# In[5]:


print(audit_df.shape)
audit_df.head(4)


# ## Interpretation of `null_threshold`
# 
# The `null_threshold` column above is also generated in [broadinstitute/cytominer_scripts/audit.R](https://github.com/broadinstitute/cytominer_scripts/blob/master/audit.R).
# 
# The metric represents a random sampling of different samples and combined, over many iterations, to represent the 95% quantile of the null distribution of pairwise correlation.

# In[6]:


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


# ## Generate Median Correlation of Replicates _Across_ Batch
# 
# We are interested in identifying the median correlation of all replicates (guides, cell lines, targets) across batches.
# 
# ### Interpretation of `median_replicate_correlation`
# 
# This value represents the median correlation across all **replicates** targeting the given gene for the given cell line **across** batches.
# 
# Therefore, the `median_replicate_correlation` is a median (guides across batches) of medians (replicates of single guides within batch) of medians (profiles of individual cells in a single well).

# In[7]:


audit_group_cols = ['Metadata_gene_name', 'Metadata_pert_name', 'Metadata_cell_line']
audit_correlation_cols = audit_group_cols + ['correlation_ci_low', 'correlation_ci_high']

median_replicate_correlation_df = (
    audit_df
    .groupby(audit_group_cols)
    .apply(get_confidence)
    .reset_index(drop=True)
    .assign(num_replicates=0)
    .groupby(audit_correlation_cols)
    .agg(
        {'correlation': 'median',
         'num_replicates': 'size'}
    )
    .reset_index()
    .rename(
        columns=
        {'correlation': 'median_replicate_correlation',
         'correlation_ci_low': 'replicate_cor_ci_low',
         'correlation_ci_high': 'replicate_cor_ci_high'}
    )
)

print(median_replicate_correlation_df.shape)
median_replicate_correlation_df.head(2)


# ## Generate Median Correlation of *Guides Targeting the Same Gene*
# 
# We are also interested in identifying the median correlation of all guides that target a single gene.
# 
# ## Interpretation of `median_target_correlation`
# 
# This value represents the median correlation across all **guides** targeting the same gene for the given cell line.
# 
# The median target correlation is calculated within batch for all profiles generated by guides targeting the same gene within each cell line.

# In[8]:


# Summarize the audits using median correlation for all guides across genes, and cell lines
audit_guide_cols = ['Metadata_gene_name', 'Metadata_cell_line']
audit_guide_correlation_cols = audit_guide_cols + ['correlation_ci_low', 'correlation_ci_high']

median_target_correlation_df = (
    median_replicate_correlation_df
    .groupby(audit_guide_cols)
    .apply(lambda x: get_confidence(x, col='median_replicate_correlation'))
    .reset_index(drop=True)
    .groupby(audit_guide_correlation_cols)
    .agg(
        {'median_replicate_correlation': 'median',
         'num_replicates': 'sum'}
    )
    .reset_index()
    .rename(
        columns=
            {'median_replicate_correlation': 'median_target_correlation',
             'num_replicates': 'num_guides',
             'correlation_ci_low': 'target_cor_ci_low',
             'correlation_ci_high': 'target_cor_ci_high'}
    )
)

print(median_target_correlation_df.shape)
median_target_correlation_df.head()


# ## Merge Audit Data Together for Plotting

# In[9]:


# Merge summarized data into audit dataframe
guide_merge_cols = ['Metadata_gene_name', 'Metadata_cell_line', 'Metadata_pert_name']
gene_merge_cols = ['Metadata_gene_name', 'Metadata_cell_line']

audit_correlation_df = (
    audit_df
    .merge(
        median_replicate_correlation_df,
        left_on=guide_merge_cols,
        right_on=guide_merge_cols,
        how='outer'
    )
    .merge(
        median_target_correlation_df,
        left_on=gene_merge_cols,
        right_on=gene_merge_cols,
        how='outer'
    )
)

# Determine how many CRISPRs are targeting the gene
audit_correlation_df = (
    audit_correlation_df
    .assign(
        guide_index=(
            audit_df
            .groupby(guide_merge_cols)
            .cumcount()+1
        )
        .tolist()
    )
    .sort_values(by='median_target_correlation',
                 ascending=False)
    .reset_index(drop=True)
)

audit_correlation_df.guide_index = audit_correlation_df.guide_index.astype('str')
audit_correlation_df.head()


# In[10]:


# Get audit_df ready to plot
gene_name_cat = pd.CategoricalDtype(categories=audit_correlation_df.Metadata_gene_name.unique().tolist())
audit_correlation_df = audit_correlation_df.assign(
    gene_plot_id=audit_correlation_df.Metadata_gene_name.astype(str).astype(gene_name_cat)
)

# Merge with the summary to add in cell line specific lines
audit_correlation_df = (
    audit_correlation_df
    .merge(summary_audit_df,
           left_on='Metadata_cell_line',
           right_on='Metadata_cell_line')
)

audit_correlation_df.head()


# ## Visualize Data
# 
# 1. Visualize the correlation of replicates within batches for each gene and cell line.
# 2. Visualize the distribution of replicates of CRISPR experiments across batches.
# 3. Visualize the correlation of CRISPR guides targeting the same gene.

# ### 1. Replicate Correlation

# In[11]:


gg.options.figure_size=(8, 15)

replicate_correlation_gg = (
    gg.ggplot(audit_correlation_df,
              gg.aes(x="correlation",
                     y="gene_plot_id")) + \
        gg.geom_jitter(size=0.5,
                       width=0,
                       height=0.4) + \
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

replicate_correlation_gg


# In[12]:


file = os.path.join("figures", "replicate-correlation")
for extension in ['.png', '.pdf']:
    gg.ggsave(replicate_correlation_gg,
              filename='{}{}'.format(file, extension),
              dpi=300,
              height=15,
              width=8,
              units='in')


# ### 2. Gene Correlation

# In[13]:


gg.options.figure_size=(6.4, 4.8)

# Make sure to drop duplicates of redundant gene, perturbation, and cell line columns
# Not removing replicates will put more weight on genes with more measurements

cor_density_gg = (
    gg.ggplot(audit_correlation_df.drop_duplicates(["Metadata_cell_line",
                                                    "Metadata_gene_name"]),
              gg.aes(x="median_target_correlation")) + \
        gg.geom_density(gg.aes(fill="Metadata_cell_line"),
                        alpha=0.4) + \
        gg.geom_vline(gg.aes(xintercept="null_threshold",
                             color="Metadata_cell_line"),
                      linetype='dashed') + \
        gg.geom_rug(gg.aes(color="Metadata_cell_line"),
                    show_legend={'color': False}) + \
        gg.theme_bw() + \
        gg.xlim([-0.5, 1]) + \
        gg.xlab("Median Correlation of All Guides Across Genes") + \
        gg.ylab("Density") + \
        gg.scale_fill_manual(name="Cell Line",
                             values=["#1b9e77", "#d95f02", "#7570b3"]) + \
        gg.scale_color_manual(name="Cell Line",
                              values=["#1b9e77", "#d95f02", "#7570b3"])
)

cor_density_gg


# In[14]:


file = os.path.join("figures", "median-guide-correlation-density")
for extension in ['.png', '.pdf']:
    gg.ggsave(cor_density_gg,
              filename='{}{}'.format(file, extension),
              dpi=300,
              height=4.8,
              width=6.4,
              units='in')


# ### 3. CRISPR Guide Correlation

# In[15]:


guide_merge_cols = ['Metadata_cell_line', 'Metadata_gene_name']

audit_guide_plot_df = (
    audit_guide_df.merge(
        median_target_correlation_df,
        left_on=guide_merge_cols,
        right_on=guide_merge_cols,
        how='outer'
    )
)

# Setup plotting logic for scatter plot
audit_guide_plot_df = audit_guide_plot_df.assign(num_guides_plot=audit_guide_plot_df.num_guides / 2)
audit_guide_plot_df.loc[audit_guide_plot_df.num_guides_plot.isna(), 'num_guides_plot'] = 1
audit_guide_plot_df.num_guides_plot = audit_guide_plot_df.num_guides_plot.astype(int)
audit_guide_plot_df.loc[audit_guide_plot_df.num_guides_plot > 5, 'num_guides_plot'] = ">5"
audit_guide_plot_df.num_guides_plot = audit_guide_plot_df.num_guides_plot.astype(str)

audit_guide_plot_df.loc[audit_guide_plot_df.num_guides_plot == '1', 'median_target_correlation'] =  audit_guide_plot_df.loc[audit_guide_plot_df.num_guides_plot == '1', 'correlation']

file = os.path.join("data", "audit_guide_replicate_correlation.tsv")
audit_guide_df.to_csv(file, sep='\t', index=False)

audit_guide_plot_df.head()


# In[16]:


get_ipython().run_cell_magic('R', '-i audit_guide_plot_df -h 3.5 -w 10 --units in -r 300', '\nsuppressPackageStartupMessages(library(dplyr))\nsuppressPackageStartupMessages(library(ggplot2))\nsuppressPackageStartupMessages(library(ggrepel))\n\naxis_title_size <- 10\naxis_text_size <- 9\nstrip_text_size <- 9\nggrepel_label_size <- 1.9\ntitle_text_size <- 10\n\ntext_color_logic <- audit_guide_plot_df$Metadata_gene_name %in% c("LacZ", "Luc", "Chr2")\ncontrol_text_color <- ifelse(text_color_logic, "red", "black")\n\naudit_guide_plot_df$num_guides_plot <-\n    factor(audit_guide_plot_df$num_guides_plot, levels = c("1", "2", "3", "4", ">5"))\n\nguide_correlation_gg <-\n    ggplot(audit_guide_plot_df,\n           aes(x=correlation,\n               y=median_target_correlation)) +\n    geom_point(aes(color=Metadata_cell_line,\n                   size=num_guides_plot),\n               alpha=0.4) +\n    geom_text_repel(arrow = arrow(length = unit(0.01, "npc")),\n                    size = ggrepel_label_size,\n                    segment.size = 0.1,\n                    segment.alpha = 0.8,\n                    force = 20,\n                    color = control_text_color,\n                    aes(x = correlation,\n                        y = median_target_correlation,\n                        label = Metadata_gene_name)) +\n    xlab("Correlation of CRISPR Guides Targeting the same Gene") +\n    ylab("Median Correlation of CRISPR Replicates") +\n    facet_grid(~Metadata_cell_line) +\n    scale_color_manual(name="Cell Line",\n                       values=c("A549" = "#1b9e77",\n                                "ES2" = "#d95f02",\n                                "HCC44" = "#7570b3")) +\n    scale_size_manual(name="Number of Guides",\n                      labels=c("1" = "1",\n                               "2" = "2",\n                               "3" = "3",\n                               "4" = "4",\n                               ">5" = ">5"),\n                      values=c("1" = 1,\n                               "2" = 2,\n                               "3" = 3,\n                               "4" = 4,\n                               ">5" = 5)) +\n    xlim(c(-0.5, 1)) +\n    theme_bw() +\n    theme(axis.text = element_text(size = axis_text_size),\n          axis.title = element_text(size = axis_title_size),\n          strip.text = element_text(size = strip_text_size),\n          strip.background = element_rect(colour = "black",\n                                          fill = "#fdfff4")) +\n    guides(color = guide_legend(order = 1))\n\nfile_base <- file.path("figures", "guide_correlation")\nfor (extension in c(\'.png\', \'.pdf\')) {\n    ggsave(guide_correlation_gg,\n           filename = paste0(file_base, extension),\n           dpi = 300,\n           height = 3.5,\n           width = 10)\n}\n\nguide_correlation_gg')

