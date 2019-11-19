#!/usr/bin/env python
# coding: utf-8

# ## Apply all Cell-Health Models to Training and Testing Sets
# 
# **Gregory Way, 2019**

# In[1]:


import os
import pandas as pd
from joblib import load

from scripts.ml_utils import load_train_test, load_models


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


def apply_model(model, feature, train_x, test_x):
    """
    Apply model to training and testing matrix
    """
    pred_train_df = (
        pd.DataFrame(model.predict(train_x), columns=["score"])
        .assign(profiles=train_x.index,
                Metadata_data_type="train",
                model=feature)
    )
    pred_test_df = (
        pd.DataFrame(model.predict(test_x), columns=["score"])
        .assign(profiles=test_x.index,
                Metadata_data_type="test",
                model=feature)
    )

    pred_df = pd.concat([pred_train_df, pred_test_df]).reset_index(drop=True)
    return pred_df

def sample_squared_error(scores, y):
    """
    Calculate the squared error per sample depending on model scores
    """
    metadata_cols = [x for x in scores.columns if x.startswith("Metadata_")]
    scores_values = scores.drop(metadata_cols, axis="columns")
    
    all_squared_error = {}
    for cell_health_feature in scores_values.columns:
        y_subset_df = y.loc[:, cell_health_feature].dropna().T
        scores_subset = scores_values.loc[:, cell_health_feature].reindex(y_subset_df.index).T

        squared_error = (y_subset_df - scores_subset) ** 2
        all_squared_error[cell_health_feature] = squared_error
    
    return pd.DataFrame(all_squared_error).reindex(scores.index)


# ## 1) Load Models and Model Coefficients
# 
# For real data and shuffled model data.

# In[4]:


model_dict, model_coef = load_models()
shuffle_model_dict, shuffle_model_coef = load_models(shuffle=True)


# In[5]:


# Load Metadata Mapping File
data_dir = os.path.join("..", "1.generate-profiles", "data")
file = os.path.join(data_dir, "profile_id_metadata_mapping.tsv")
metadata_df = pd.read_csv(file, sep='\t')

metadata_df.head()


# ## 2) Load Training and Testing Data

# In[6]:


x_train_df, x_test_df, y_train_df, y_test_df = load_train_test(drop_metadata=True)


# ## 3) Output Model Coefficients

# In[7]:


# Extract all model coefficients and output to file
coef_df = pd.DataFrame(model_coef)
coef_df.index = x_test_df.columns
coef_df.index.name = "features"

file = os.path.join("results", "all_model_coefficients.tsv")
coef_df.to_csv(file, sep='\t', index=True)

print(coef_df.shape)
coef_df.head(2)


# In[8]:


# Extract all model coefficients and output to file
shuffle_coef_df = pd.DataFrame(shuffle_model_coef)
shuffle_coef_df.index = x_test_df.columns
shuffle_coef_df.index.name = "features"

file = os.path.join("results", "all_model_coefficients_shuffled.tsv")
shuffle_coef_df.to_csv(file, sep='\t', index=True)

print(shuffle_coef_df.shape)
shuffle_coef_df.head(2)


# ## 4) Apply all models
# 
# For real and shuffled data.

# In[9]:


all_scores = []
all_shuffle_scores = []
for cell_health_feature in model_dict.keys():
    # Apply Real Model Classifiers
    model_clf = model_dict[cell_health_feature]
    pred_df = apply_model(model=model_clf,
                          feature=cell_health_feature,
                          train_x=x_train_df,
                          test_x=x_test_df)
    all_scores.append(pred_df)
    
    # Apply Shuffled Model Classifiers
    shuffle_model_clf = shuffle_model_dict[cell_health_feature]
    shuffle_pred_df = apply_model(model=shuffle_model_clf,
                                  feature=cell_health_feature,
                                  train_x=x_train_df,
                                  test_x=x_test_df)
    all_shuffle_scores.append(shuffle_pred_df)


# ## 5) Concatenate scores with Metadata

# In[10]:


# Concatenate real data scores
all_scores = (
    pd.concat(all_scores)
    .reset_index(drop=True)
    .pivot_table(index=["profiles", "Metadata_data_type"],
                 columns="model",
                 values="score")
    .reset_index()
)

all_scores = (
    metadata_df.merge(all_scores,
                      left_on="Metadata_profile_id",
                      right_on="profiles")
    .drop("profiles", axis="columns")
)

all_scores.index = all_scores.Metadata_profile_id
all_scores = all_scores.drop("Metadata_profile_id", axis="columns")

# Output file
file = os.path.join("results", "all_model_predictions.tsv")
all_scores.to_csv(file, sep='\t', index=True)

print(all_scores.shape)
all_scores.head(2)


# In[11]:


# Concatenate shuffled data scores
all_shuffle_scores = (
    pd.concat(all_shuffle_scores)
    .reset_index(drop=True)
    .pivot_table(index=["profiles", "Metadata_data_type"],
                 columns="model",
                 values="score")
    .reset_index()
)

all_shuffle_scores = (
    metadata_df.merge(all_shuffle_scores,
                      left_on="Metadata_profile_id",
                      right_on="profiles")
    .drop("profiles", axis="columns")
)

all_shuffle_scores.index = all_shuffle_scores.Metadata_profile_id
all_shuffle_scores = all_shuffle_scores.drop("Metadata_profile_id", axis="columns")

# Output file
file = os.path.join("results", "all_model_predictions_shuffled.tsv")
all_shuffle_scores.to_csv(file, sep='\t', index=True)

print(all_shuffle_scores.shape)
all_shuffle_scores.head(2)


# ## 6) Calculate the Squared Error of Individual Samples
# 
# For real and shuffled data

# In[12]:


y_df = pd.concat([y_train_df, y_test_df]).reindex(all_scores.index)

print(y_df.shape)
y_df.head(2)


# In[13]:


all_score_error = sample_squared_error(scores=all_scores, y=y_df)

all_score_error = (
    metadata_df.merge(all_score_error,
                      left_on="Metadata_profile_id",
                      right_index=True)
)

# Output file
file = os.path.join("results", "all_model_sample_squared_error.tsv")
all_score_error.to_csv(file, sep='\t', index=False)

print(all_score_error.shape)
all_score_error.head(2)


# In[14]:


all_shuffle_score_error = sample_squared_error(scores=all_shuffle_scores, y=y_df)

all_shuffle_score_error = (
    metadata_df.merge(all_shuffle_score_error,
                      left_on="Metadata_profile_id",
                      right_index=True)
)

# Output file
file = os.path.join("results", "all_model_sample_squared_error_shuffled.tsv")
all_shuffle_score_error.to_csv(file, sep='\t', index=False)

print(all_shuffle_score_error.shape)
all_shuffle_score_error.head()

