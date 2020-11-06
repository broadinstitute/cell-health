#!/usr/bin/env python
# coding: utf-8

# ## Predicting Cell Health Phenotypes with Specific Feature Subsets
# 
# The contribution of specific feature groups and Cell Painting stains to model performance is a major question.
# 
# Here, we address this question by systematically removing stains and feature groups and observing impact on model performance.

# In[1]:


import sys
import pathlib
import numpy as np
import pandas as pd

from sklearn.linear_model import SGDClassifier, ElasticNet
from sklearn.pipeline import Pipeline
from sklearn.exceptions import ConvergenceWarning

sys.path.insert(0, "../3.train/scripts/")
from ml_utils import CellHealthPredict, load_train_test


# In[2]:


# Set constants
consensus = "modz"
data_dir = pathlib.Path("../3.train/data")
results_dir = pathlib.Path("../3.train/results")

coef_file = pathlib.Path(f"{results_dir}/all_model_coefficients_{consensus}.tsv")

output_dir = pathlib.Path("results")
output_dir.mkdir(exist_ok=True)
output_file = pathlib.Path(f"{output_dir}/feature_group_subset_removal_robustness_results_{consensus}.tsv.gz")

# Define which features to remove per group
feature_removal_ids = {
    "feature_group": [
        "Texture",
        "Intensity",
        "Correlation",
        "AreaShape",
        "Granularity",
        "RadialDistribution",
        "Neighbors"
    ],
    "channel": [
        "AGP",
        "DNA",
        "ER",
        "Mito",
        "RNA"
    ],
    "compartment": [
        "Cells",
        "Cytoplasm",
        "Nuclei"
    ]
}


# In[3]:


# Get all features and define feature groups
feature_df = pd.read_csv(coef_file, sep="\t")

feature_info_split_df = (
    pd.concat(
        [pd.Series(x.split("_")) for x in feature_df.features],
        axis="columns"
    )
    .transpose()
)
feature_info_split_df.columns = [
    "compartment", "feature_group", "measurement", "channel", "parameter1", "parameter2"
]
feature_info_split_df = feature_info_split_df.assign(full_feature_name=feature_df.features)

print(feature_info_split_df.shape)
feature_info_split_df.head()


# In[4]:


# Set ML constants
# We will optimize each model independently
alphas = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
l1_ratios = [0.1, 0.12, 0.14, 0.16, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9]
n_folds = 5

parameters = {
    'regress__alpha': alphas,
    'regress__l1_ratio': l1_ratios
}

estimator = Pipeline(
    steps=[(
        "regress",
        ElasticNet(
            random_state=42,
            max_iter=50,
            tol=1e-3
        )
        
    )]
)

y_transform = "raw"
scoring = "r2"
decision_function = False
model_type = "Regression"


# In[5]:


# Load training/test data
x_train_df, x_test_df, y_train_df, y_test_df = (
    load_train_test(data_dir=data_dir, drop_metadata=True, consensus=consensus)
)

# Setup target variable names and assert that there are 70 of them
cell_health_targets = y_train_df.columns.tolist()
assert len(cell_health_targets) == 70

print(x_train_df.shape)
x_train_df.head(3)


# ## Systematically drop features and save performance results

# In[6]:


regression_results_list = []
for feature_category in feature_removal_ids:
    feature_removal = feature_removal_ids[feature_category]
    for feature_type in feature_removal:
        subset_features_to_drop = (
            feature_info_split_df
            .loc[
                feature_info_split_df.loc[:, feature_category] == feature_type,
                "full_feature_name"
            ]
            .tolist()
        )
        
        num_features_dropped = len(subset_features_to_drop)
        x_train_subset_df = x_train_df.drop(subset_features_to_drop, axis="columns")
        x_test_subset_df = x_test_df.drop(subset_features_to_drop, axis="columns")
        
        for cell_health_target in cell_health_targets:
            print(f"Now Training Target: {cell_health_target}")
            print(f"[Class] Sample Group Removal: {feature_category}; Dropping: {feature_type}\n")
            
            # Initialize class
            chp = CellHealthPredict(
                x_df=x_train_subset_df,
                y_df=y_train_df,
                parameters=parameters,
                estimator=estimator,
                n_folds=n_folds,
                cv_scoring=scoring,
                shuffle=False
            )

            # Fit model
            is_fit = chp.fit_cell_health_target(
                cell_health_target,
                y_transform=y_transform
            )

            # Training performance metrics
            metric_mse, metric_rtwo, y_true, y_pred = chp.get_performance(
                decision_function=decision_function,
                return_y=True,
                cell_line="all"
            )

            # Testing performance metrics
            metric_mse_test, metric_rtwo_test, y_test_true, y_test_pred = chp.get_performance(
                x_test=x_test_subset_df,
                y_test=y_test_df,
                decision_function=decision_function,
                return_y=True,
                data_fit_type="test",
                cell_line="all"
            )

            # Store results
            mse_df = (
                pd.concat([metric_mse, metric_mse_test], axis='rows')
                .reset_index(drop=True)
                .assign(
                    feature_category=feature_category,
                    feature_type_dropped=feature_type,
                    num_dropped=num_features_dropped
                )
            )
            
            rtwo_df = (
                pd.concat([metric_rtwo, metric_rtwo_test], axis='rows')
                .reset_index(drop=True)
                .assign(
                    feature_category=feature_category,
                    feature_type_dropped=feature_type,
                    num_dropped=num_features_dropped
                )
            )
            
            regression_results_list.append(mse_df)
            regression_results_list.append(rtwo_df)


# In[7]:


# Compile the full regression results
full_regression_results_df = pd.concat(regression_results_list).reset_index(drop=True)

print(full_regression_results_df.shape)
full_regression_results_df.head(3)


# In[8]:


# Save all results
full_regression_results_df.to_csv(output_file, sep="\t", index=False, compression="gzip")

