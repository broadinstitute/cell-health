#!/usr/bin/env python
# coding: utf-8

# # Compare cytominer profiles to pycytominer profiles
# 
# **Gregory Way, 2019**

# In[1]:


import os
import pandas as pd
from sklearn.metrics import mean_squared_error
import plotnine as gg


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# ## Load Profiles
# 
# ### 1) Profiles processed using cytominer

# In[3]:


plate = "SQ00014612"

cyto_file = os.path.join("profiles", "{}_normalized_variable_selected.csv".format(plate))
cyto_df = pd.read_csv(cyto_file)

print(cyto_df.shape)
cyto_df.head(2)


# ### 2) Profiles processed using pycytominer

# In[4]:


batch = "CRISPR_PILOT_B1"
pycyto_base = os.path.join("..", "..", "1.generate-profiles", "data", "profiles", batch, plate)
pycyto_file = os.path.join(pycyto_base, "{}_normalized_feature_select.csv.gz".format(plate))
pycyto_df = pd.read_csv(pycyto_file).rename({"Image_Metadata_Plate": "Metadata_Plate",
                                             "Image_Metadata_Well": "Metadata_Well"},
                                            axis='columns')

print(pycyto_df.shape)
pycyto_df.head(2)


# ## Subset CellProfiler Features

# In[5]:


meta_cols = ["Metadata_Plate", "Metadata_Well", "Metadata_gene_name", "Metadata_pert_name"]
cp_features = meta_cols + [
    x for x in cyto_df.columns if
    (x.startswith("Cells_") or x.startswith("Nuclei_") or x.startswith("Cytoplasm_"))
]

len(cp_features)


# In[6]:


cyto_df = (
    cyto_df
    .reindex(cp_features, axis="columns")
    .sort_values(by="Metadata_Well")
    .assign(Metadata_Package="cytominer")
)
cyto_df.head(2)


# In[7]:


pycyto_df = (
    pycyto_df
    .reindex(cp_features, axis="columns")
    .sort_values(by="Metadata_Well")
    .assign(Metadata_Package="pycytominer")
)
pycyto_df.head(2)


# In[8]:


meta_cols = meta_cols + ["Metadata_Package"]

df = pd.concat([cyto_df, pycyto_df]).reset_index(drop=True)

print(df.shape)
df.head(2)


# In[9]:


cp_df = df.drop(meta_cols, axis="columns")
meta_df = df.drop(cp_df.columns, axis="columns")
meta_df = meta_df.assign(meta_id = meta_df.Metadata_Well + "_" + meta_df.Metadata_Package)


# In[10]:


cp_cor_df = cp_df.transpose().corr()
cp_cor_df = cp_cor_df.where(pd.np.tril(pd.np.ones(cp_cor_df.shape), k=-1).astype(bool))

cp_cor_df.index = meta_df.Metadata_Well + "_" + meta_df.Metadata_Package
cp_cor_df.columns = meta_df.Metadata_Well + "_" + meta_df.Metadata_Package
cp_cor_df = cp_cor_df.stack().reset_index().rename({"level_0": "a", "level_1": "b", 0: "correlation"},
                                                   axis='columns')

print(cp_cor_df.shape)
cp_cor_df.head(2)


# In[11]:


full_df = (
    cp_cor_df
    .merge(meta_df,
           left_on="a",
           right_on="meta_id")
    .merge(meta_df,
           left_on="b",
           right_on="meta_id",
           suffixes=["", "_b"])
    .drop(["Metadata_Plate", "Metadata_Plate_b",
           "meta_id", "meta_id_b"], axis="columns")
)

full_df = full_df.assign(
    across_package=full_df.Metadata_Well == full_df.Metadata_Well_b,
    across_package_gene=(
        (full_df.Metadata_Package != full_df.Metadata_Package_b) &
        (full_df.Metadata_gene_name == full_df.Metadata_gene_name_b) &
        (full_df.Metadata_Well != full_df.Metadata_Well_b)
    ),
    within_package_gene=(
        (full_df.Metadata_Package == full_df.Metadata_Package_b) &
        (full_df.Metadata_gene_name == full_df.Metadata_gene_name_b)
    )
)

full_df.loc[full_df.across_package_gene == True, "across_package"] = "Same Gene Across Package"
full_df.loc[full_df.within_package_gene == True, "across_package"] = "Same Gene Within Package"
full_df.loc[full_df.across_package == True, "across_package"] = "Same Well Across Package"
full_df.loc[full_df.across_package == False, "across_package"] = "Different Well"

print(full_df.shape)
full_df.head(2)


# ## Visualize Correlation Distribution

# In[12]:


cor_gg = gg.ggplot(full_df, gg.aes(x="correlation")) +     gg.geom_density(gg.aes(fill="across_package"), alpha=0.8) +     gg.theme_bw() +     gg.scale_fill_discrete(name="") +     gg.facet_wrap("across_package", scales='free_y', nrow=4)

fig_file = os.path.join("compare_cytominer_pycytominer_pearson.png")
cor_gg.save(filename=fig_file, dpi=300, width=5, height=5)

cor_gg


# ## Calculate pairwise MSE

# In[13]:


all_mse_list = []
n = df.shape[0]
for i in range(0, n):
    for j in range(0, n):
        a = meta_df.iloc[i, :].rename('{}_a'.format)
        b = meta_df.iloc[j, :].rename('{}_b'.format)
        mse_result = mean_squared_error(cp_df.iloc[i, :].fillna(0), cp_df.iloc[j, :].fillna(0))
        all_mse_list.append(pd.concat([a, b, pd.Series(mse_result, index=["mse"])]))


# In[14]:


all_mse_df = pd.concat(all_mse_list, axis="columns").transpose()

print(all_mse_df.shape)
all_mse_df.head()


# In[15]:


all_mse_df = all_mse_df.assign(
    across_package=all_mse_df.Metadata_Well_a == all_mse_df.Metadata_Well_b,
    across_package_gene=(
        (all_mse_df.Metadata_Package_a != all_mse_df.Metadata_Package_b) &
        (all_mse_df.Metadata_gene_name_a == all_mse_df.Metadata_gene_name_b) &
        (all_mse_df.Metadata_Well_a != all_mse_df.Metadata_Well_b)
    ),
    within_package_gene=(
        (all_mse_df.Metadata_Package_a == all_mse_df.Metadata_Package_b) &
        (all_mse_df.Metadata_gene_name_a == all_mse_df.Metadata_gene_name_b)
    )
)

all_mse_df.loc[all_mse_df.across_package_gene == True, "across_package"] = "Same Gene Across Package"
all_mse_df.loc[all_mse_df.within_package_gene == True, "across_package"] = "Same Gene Within Package"
all_mse_df.loc[all_mse_df.across_package == True, "across_package"] = "Same Well Across Package"
all_mse_df.loc[all_mse_df.across_package == False, "across_package"] = "Different Well"

all_mse_df.mse = all_mse_df.mse.astype(float)
all_mse_df.across_package = all_mse_df.across_package.astype(str)

all_mse_df = all_mse_df.drop_duplicates(subset=["mse"])

print(all_mse_df.shape)
all_mse_df.head()


# In[16]:


mse_gg = gg.ggplot(all_mse_df, gg.aes(x="mse")) +     gg.geom_density(gg.aes(fill="across_package"), alpha=0.2) +     gg.theme_bw() +     gg.scale_fill_discrete(name="") +     gg.xlim([0, 5])

fig_file = os.path.join("compare_cytominer_pycytominer_mse.png")
mse_gg.save(filename=fig_file, dpi=300, width=5, height=5)

mse_gg

