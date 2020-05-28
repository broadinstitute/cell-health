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
import pandas as pd

import matplotlib.pyplot as plt
import plotnine as gg

from pycytominer import audit

from scripts.audit_utils import get_confidence


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[3]:


profile_file = os.path.join(
    "..", "1.generate-profiles", "data", "processed", "cell_health_profiles_merged.tsv.gz"
)

profile_df = pd.read_csv(profile_file, sep='\t')

print(profile_df.shape)
profile_df.head()


# ## Perform Audit of Genes and Guides 

# In[4]:


audit_gene_groups = ["Metadata_cell_line", "Metadata_gene_name"]
audit_gene_df = audit(
    profile_df,
    audit_groups=audit_gene_groups,
    iterations=100
)

print(audit_gene_df.shape)
audit_gene_df.head()


# In[5]:


audit_guide_groups = ["Metadata_cell_line", "Metadata_gene_name", "Metadata_pert_name"]
audit_guide_df = audit(
    profile_df,
    audit_groups=audit_guide_groups,
    iterations=100
)

print(audit_guide_df.shape)
audit_guide_df.head()


# ## Get Median Same Gene Guide Correlation

# In[6]:


same_gene_groupby_cols = audit_gene_groups + ["replicate_type"]

median_cor_across_same_gene_guides_df = (
    audit_guide_df
    .groupby(
        same_gene_groupby_cols
    )["correlation"]
    .median()
    .reset_index()
    .rename(
        {
            "correlation": "median_same_gene_guide_correlation"
        },
        axis="columns"
    )
)

print(median_cor_across_same_gene_guides_df.shape)
median_cor_across_same_gene_guides_df.head()


# ## Count Number of Unique Guides per Gene

# In[7]:


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

# In[8]:


summary_corr_df = (
    audit_gene_df
    .merge(
        median_cor_across_same_gene_guides_df,
        on=same_gene_groupby_cols
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
            "replicate": "Replicate",
            "non_replicate": "Non-Replicate"
        }
    )
)

print(summary_corr_df.shape)
summary_corr_df.head()


# ## Generate Summary Figures

# In[9]:


get_ipython().run_cell_magic('R', '-i summary_corr_df -h 3.5 -w 10 --units in -r 300', '\nsuppressPackageStartupMessages(library(dplyr))\nsuppressPackageStartupMessages(library(ggplot2))\nsuppressPackageStartupMessages(library(ggrepel))\n\nsource(file.path("..", "3.train", "scripts", "assay_themes.R"))\n\naxis_title_size <- 10\naxis_text_size <- 9\nstrip_text_size <- 9\nggrepel_label_size <- 1.9\ntitle_text_size <- 10\n\nsummary_corr_df$num_guides_plot <-\n    factor(summary_corr_df$num_guides_plot, levels = c("1", "2", "3", ">4"))\n\naudit_guide_plot_df <- summary_corr_df %>% dplyr::filter(replicate_type == "Replicate")\n\ntext_color_logic <- audit_guide_plot_df$Metadata_gene_name %in% c("LacZ", "Luc", "Chr2")\ncontrol_text_color <- ifelse(text_color_logic, "red", "black")\n\nguide_correlation_gg <-\n    ggplot(audit_guide_plot_df,\n           aes(x=correlation,\n               y=median_same_gene_guide_correlation)) +\n    geom_point(aes(color=Metadata_cell_line,\n                   size=num_guides_plot),\n               alpha=0.4) +\n    geom_text_repel(arrow = arrow(length = unit(0.01, "npc")),\n                    size = ggrepel_label_size,\n                    segment.size = 0.1,\n                    segment.alpha = 0.8,\n                    force = 20,\n                    color = control_text_color,\n                    aes(x = correlation,\n                        y = median_same_gene_guide_correlation,\n                        label = Metadata_gene_name)) +\n    xlab("Correlation of CRISPR Guides Targeting the same Gene") +\n    ylab("Median Correlation of CRISPR Replicates") +\n    facet_grid(~Metadata_cell_line) +\n    scale_color_manual(name="Cell Line",\n                       values=cell_line_colors,\n                       labels=cell_line_labels) +\n    scale_size_manual(name="Number of Guides",\n                      labels=c("1" = "1",\n                               "2" = "2",\n                               "3" = "3",\n                               ">4" = ">4"),\n                      values=c("1" = 1,\n                               "2" = 2,\n                               "3" = 3,\n                               ">4" = 4)) +\n    xlim(c(-0.5, 1)) +\n    theme_bw() +\n    theme(\n        axis.text = element_text(size = axis_text_size),\n        axis.title = element_text(size = axis_title_size),\n        strip.text = element_text(size = strip_text_size),\n        strip.background = element_rect(colour = "black", fill = "#fdfff4")\n    ) +\n    guides(color = guide_legend(order = 1))\n\nfile_base <- file.path("figures", "guide_correlation")\nfor (extension in c(\'.png\', \'.pdf\')) {\n    ggsave(guide_correlation_gg,\n           filename = paste0(file_base, extension),\n           dpi = 300,\n           height = 3.5,\n           width = 10)\n}\n\nguide_correlation_gg')


# In[10]:


gg.options.figure_size=(6.4, 4.8)

# Make sure to drop duplicates of redundant gene, perturbation, and cell line columns
# Not removing replicates will put more weight on genes with more measurements

cor_density_gg = (
    gg.ggplot(
        summary_corr_df.drop_duplicates(
            ["Metadata_cell_line", "Metadata_gene_name", "replicate_type"]
        ),
        gg.aes(x="median_same_gene_guide_correlation")) + \
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

