#!/usr/bin/env python
# coding: utf-8

# # Analyzing Audits
# 
# **Predicting Cell Health**
# 
# **Gregory Way, 2020**
# 
# Detecting replicate correlation across replicates and CRISPR guides targeting the same gene.

# In[1]:


import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import plotnine as gg

from pycytominer import audit

from scripts.audit_utils import get_confidence


# In[2]:


np.random.seed(123)


# In[3]:


get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[4]:


profile_file = os.path.join(
    "..", "1.generate-profiles", "data", "processed", "cell_health_profiles_merged.tsv.gz"
)

profile_df = pd.read_csv(profile_file, sep='\t')

print(profile_df.shape)
profile_df.head()


# In[5]:


# How many replicates per perturbation?
(
    profile_df
    .groupby(["Metadata_cell_line","Metadata_gene_name", "Metadata_pert_name"])
    ["Metadata_Plate"]
    .count()
    .value_counts()
)


# ## Perform Audit of Genes and Guides 

# In[6]:


audit_gene_groups = ["Metadata_cell_line", "Metadata_gene_name"]
audit_guide_groups = ["Metadata_cell_line", "Metadata_gene_name", "Metadata_pert_name"]

# Perform a full audit and collect cell line, gene, and guide metadata
audit_gene_df = audit(
    profile_df,
    audit_groups=audit_guide_groups,
    audit_resolution="full"
)

print(audit_gene_df.shape)
audit_gene_df.head()


# In[7]:


# Determine gene replicates
pair_a_df = audit_gene_df.loc[:, ["{}_pair_a".format(x) for x in audit_gene_groups]]
pair_a_df.columns = audit_gene_groups
pair_b_df = audit_gene_df.loc[:, ["{}_pair_b".format(x) for x in audit_gene_groups]]
pair_b_df.columns = audit_gene_groups

audit_gene_df = audit_gene_df.assign(
    replicate_gene_info=(pair_a_df == pair_b_df).all(axis="columns")
)

# Determine guide replicates
pair_a_df = audit_gene_df.loc[:, ["{}_pair_a".format(x) for x in audit_guide_groups]]
pair_a_df.columns = audit_guide_groups
pair_b_df = audit_gene_df.loc[:, ["{}_pair_b".format(x) for x in audit_guide_groups]]
pair_b_df.columns = audit_guide_groups

audit_gene_df = audit_gene_df.assign(
    replicate_guide_info=(pair_a_df == pair_b_df).all(axis="columns")
)

print(audit_gene_df.shape)
audit_gene_df.head()


# In[8]:


pd.crosstab(
    audit_gene_df.replicate_gene_info,
    audit_gene_df.replicate_guide_info
)


# In[9]:


gene_df = (
    audit_gene_df
    .query("not replicate_guide_info")
    .groupby(["{}_pair_a".format(x) for x in audit_gene_groups] + ["replicate_gene_info"])
    ["pairwise_correlation"]
    .median()
    .reset_index()
    .rename(
        {
            "Metadata_cell_line_pair_a": "Metadata_cell_line",
            "Metadata_gene_name_pair_a": "Metadata_gene_name",
            "replicate_gene_info": "replicate_type",
            "pairwise_correlation": "correlation"
        },
        axis="columns"
    )
)

print(gene_df.shape)
gene_df.head()


# In[10]:


guide_df = (
    audit_gene_df
    .groupby(["{}_pair_a".format(x) for x in audit_gene_groups] + ["replicate_guide_info"])
    ["pairwise_correlation"]
    .median()
    .reset_index()
    .rename(
        {
            "Metadata_cell_line_pair_a": "Metadata_cell_line",
            "Metadata_gene_name_pair_a": "Metadata_gene_name",
            "replicate_guide_info": "replicate_type",
            "pairwise_correlation": "correlation"
        },
        axis="columns"
    )
)

print(guide_df.shape)
guide_df.head()


# In[11]:


gene_non_replicate_quantile = (
    pd.DataFrame(
        gene_df
        .query("not replicate_type")
        .groupby("Metadata_cell_line")["correlation"]
        .quantile(0.95)
    )
    .reset_index()
    .assign(random="95% Non-Replicate")
    .rename({"correlation": "correlation_null"}, axis = "columns")
)

gene_non_replicate_quantile


# In[12]:


guide_non_replicate_quantile = (
    pd.DataFrame(
        guide_df
        .query("not replicate_type")
        .groupby("Metadata_cell_line")["correlation"]
        .quantile(0.95)
    )
    .reset_index()
    .assign(random="95% Non-Replicate")
    .rename({"correlation": "correlation_null"}, axis = "columns")
)

guide_non_replicate_quantile


# In[13]:


# Figure out the percentage of "strong" phenotypes
gene_strong = (
    gene_df
    .merge(gene_non_replicate_quantile, on="Metadata_cell_line", how="left")
    .query("replicate_type")
)

gene_strong = gene_strong.assign(gene_strong=gene_strong.correlation > gene_strong.correlation_null)

gene_strong = pd.DataFrame(
    gene_strong.groupby("Metadata_cell_line")["gene_strong"].sum() /
    gene_strong.groupby("Metadata_cell_line")["gene_strong"].count()
).reset_index()

# Strong phenotypes per guide profile
guide_strong = (
    guide_df
    .merge(guide_non_replicate_quantile, on="Metadata_cell_line", how="left")
    .query("replicate_type")
)

guide_strong = guide_strong.assign(guide_strong=guide_strong.correlation > guide_strong.correlation_null)

guide_strong = pd.DataFrame(
    guide_strong.groupby("Metadata_cell_line")["guide_strong"].sum() /
    guide_strong.groupby("Metadata_cell_line")["guide_strong"].count()
).reset_index()

# Combine Gene and Guide Strength to Plot Text
profile_strong = gene_strong.merge(guide_strong, on="Metadata_cell_line")
profile_strong = profile_strong.assign(
    guide_strong_text=(
        "Guide Strong: " +
        (profile_strong.guide_strong * 100).round(2).astype(str) + "%"
    ),
    gene_strong_text=(
        "Gene Strong: " +
        (profile_strong.gene_strong * 100).round(2).astype(str) + "%"
    )
)

profile_strong = profile_strong.assign(
    strength_text=profile_strong.gene_strong_text + "\n" + profile_strong.guide_strong_text
)
profile_strong


# In[14]:


profile_strong.mean()


# ## Count Number of Unique Guides per Gene

# In[15]:


guide_count_df = (
    profile_df
    .drop_duplicates(subset=audit_guide_groups)
    .groupby(audit_gene_groups)
    ["Metadata_Plate"]
    .count()
    .reset_index()
    .rename(
        {
            "Metadata_Plate": "num_unique_guides"
        },
        axis="columns"
    )
)

print(guide_count_df.shape)
guide_count_df.head()


# ## Process Data for Plotting

# In[16]:


summary_corr_df = (
    gene_df
    .merge(
        guide_df,
        on=audit_gene_groups + ["replicate_type"],
        suffixes=["_gene", "_guide"]
    )
    .merge(
        guide_count_df,
        on=audit_gene_groups
    )
)

summary_corr_df = summary_corr_df.assign(num_guides_plot=summary_corr_df.num_unique_guides)
summary_corr_df.loc[summary_corr_df.num_guides_plot.isna(), 'num_guides_plot'] = 1
summary_corr_df.num_guides_plot = summary_corr_df.num_guides_plot.astype(int)
summary_corr_df.loc[summary_corr_df.num_guides_plot > 4, 'num_guides_plot'] = ">4"
summary_corr_df.num_guides_plot = summary_corr_df.num_guides_plot.astype(str)

summary_corr_df.replicate_type = (
    summary_corr_df
    .replicate_type
    .replace(
        {
            True: "Replicate",
            False: "Non-Replicate"
        }
    )
)

print(summary_corr_df.shape)
summary_corr_df.head()


# ## Generate Summary Figures

# In[17]:


get_ipython().run_cell_magic('R', '-i summary_corr_df -i guide_non_replicate_quantile -i gene_non_replicate_quantile -i profile_strong -h 3.5 -w 10 --units in -r 300', '\nsuppressPackageStartupMessages(library(dplyr))\nsuppressPackageStartupMessages(library(ggplot2))\nsuppressPackageStartupMessages(library(ggrepel))\n\nsource(file.path("..", "3.train", "scripts", "assay_themes.R"))\n\naxis_title_size <- 10\naxis_text_size <- 9\nstrip_text_size <- 9\nggrepel_label_size <- 1.9\ntitle_text_size <- 10\n\nsummary_corr_df$num_guides_plot <-\n    factor(summary_corr_df$num_guides_plot, levels = c("1", "2", "3", ">4"))\n\naudit_guide_plot_df <- summary_corr_df %>% dplyr::filter(replicate_type == "Replicate")\n\ntext_color_logic <- audit_guide_plot_df$Metadata_gene_name %in% c("LacZ", "Luc", "Chr2")\ncontrol_text_color <- ifelse(text_color_logic, "red", "black")\n\nguide_correlation_gg <-\n    ggplot(audit_guide_plot_df,\n           aes(x = correlation_gene,\n               y = correlation_guide)) +\n    geom_abline(yintercept = 0, slope = 1, linetype = "dashed", color = "blue") +\n    geom_vline(data = gene_non_replicate_quantile,\n               aes(xintercept = correlation_null),\n               color = "red",\n               linetype = "dashed") +\n    geom_hline(data = guide_non_replicate_quantile,\n               aes(yintercept = correlation_null, linetype = random),\n               color = "red") +\n    geom_point(aes(color=Metadata_cell_line,\n                   size=num_guides_plot),\n               alpha = 0.4) +\n    geom_text(data = profile_strong,\n              x = -0.16,\n              y = 0.875,\n              size = 3,\n              aes(label = strength_text)) +\n    geom_text_repel(arrow = arrow(length = unit(0.01, "npc")),\n                    size = ggrepel_label_size,\n                    segment.size = 0.1,\n                    segment.alpha = 0.8,\n                    force = 20,\n                    color = control_text_color,\n                    aes(x = correlation_gene,\n                        y = correlation_guide,\n                        label = Metadata_gene_name)) +\n    xlab("Median Correlation of CRISPR Guides Targeting the same Gene (Remove Same Guide)") +\n    ylab("Median Correlation of CRISPR Replicates") +\n    facet_grid(~Metadata_cell_line) +\n    scale_color_manual(name = "Cell Line",\n                       values = cell_line_colors,\n                       labels = cell_line_labels) +\n    scale_size_manual(name = "Number of Guides",\n                      labels=c("1" = "1",\n                               "2" = "2",\n                               "3" = "3",\n                               ">4" = ">4"),\n                      values=c("1" = 1,\n                               "2" = 2,\n                               "3" = 3,\n                               ">4" = 4)) +\n    scale_linetype_manual(name = "Null Distribution", values = "dashed") +\n    xlim(c(-0.5, 1)) +\n    theme_bw() +\n    theme(\n        axis.text = element_text(size = axis_text_size),\n        axis.title = element_text(size = axis_title_size),\n        strip.text = element_text(size = strip_text_size),\n        strip.background = element_rect(colour = "black", fill = "#fdfff4")\n    ) +\n    guides(color = guide_legend(order = 1),\n           size = guide_legend(order = 2),\n           linetype = guide_legend(order = 3))\n\nfile_base <- file.path("figures", "guide_correlation")\nfor (extension in c(\'.png\', \'.pdf\')) {\n    ggsave(guide_correlation_gg,\n           filename = paste0(file_base, extension),\n           dpi = 300,\n           height = 3.5,\n           width = 10)\n}\n\nguide_correlation_gg')


# In[18]:


gg.options.figure_size=(6.4, 4.8)

# Make sure to drop duplicates of redundant gene, perturbation, and cell line columns
# Not removing replicates will put more weight on genes with more measurements

cor_density_gg = (
    gg.ggplot(
        summary_corr_df.drop_duplicates(
            ["Metadata_cell_line", "Metadata_gene_name", "replicate_type"]
        ),
        gg.aes(x="correlation_guide")) + \
        gg.geom_density(gg.aes(fill="Metadata_cell_line"),
                        alpha=0.4) + \
        gg.geom_rug(gg.aes(color="Metadata_cell_line"),
                    show_legend={'color': False}) + \
        gg.theme_bw() + \
    gg.theme(
            subplots_adjust={"wspace": 0.2},
            axis_text=gg.element_text(size=7),
            axis_title=gg.element_text(size=9),
            strip_text=gg.element_text(size=6, color="black"),
            strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
        ) + \
        gg.xlim([-0.5, 1]) + \
        gg.xlab("Median Correlation of All Guides Across Genes") + \
        gg.ylab("Density") + \
        gg.facet_wrap("~replicate_type", nrow=2, scales="free") + \
        gg.scale_fill_manual(name="Cell Line",
                             values=["#1b9e77", "#d95f02", "#7570b3"]) + \
        gg.scale_color_manual(name="Cell Line",
                              values=["#1b9e77", "#d95f02", "#7570b3"])
)

file = os.path.join("figures", "median-guide-correlation-density")
for extension in ['.png', '.pdf']:
    gg.ggsave(
        cor_density_gg,
        filename='{}{}'.format(file, extension),
        dpi=500,
        height=2,
        width=3,
        units='in'
    )

cor_density_gg

