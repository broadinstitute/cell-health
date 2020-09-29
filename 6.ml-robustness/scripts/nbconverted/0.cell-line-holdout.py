#!/usr/bin/env python
# coding: utf-8

# ## Performing a Cell Line holdout analysis
# 
# We used three cell lines to train our cell health models: A549, ES2, and HCC44.
# 
# Perform the following analysis:
# 
# 1. Select 2 of 3 cell lines
# 2. Train a model with all samples from the 2 cell lines
# 3. Apply the trained model to the remaining 1 cell line not used for training
# 4. Extract out R2 metrics
# 5. Perform this procedure for all 70 cell health models
# 6. Perform this procedure holding out each cell line once.
# 
# ### Data splits
# 
# 1. Train: A549, ES2; Test: HCC44
# 2. Train: A549, HCC44; Test: ES2
# 3. Train: ES2, HCC44; Test: A549
# 
# ### Knowledge gain
# 
# This analysis will tell us the extent to which certain models can be predicted in cell lines entirely unseen by the training procedure, and how well we could expect these models to generalize to new cell lines.

# In[1]:


import sys
import pathlib
import pandas as pd

from sklearn.linear_model import SGDClassifier, ElasticNet
from sklearn.pipeline import Pipeline
from sklearn.exceptions import ConvergenceWarning

sys.path.insert(0, "../3.train/scripts/")
from ml_utils import CellHealthPredict


# In[2]:


def process_cell_line_train_test(df, cell_line, meta_cols_to_drop, test_set=False):
    
    if test_set:
        output_df = df.query("Metadata_cell_line == @cell_line")
    else:
        output_df = df.query("Metadata_cell_line != @cell_line")
        
    output_df = (
        output_df
        .drop(meta_cols_to_drop, axis="columns")
        .reset_index(drop=True)
    )
    
    return output_df


# In[3]:


# Set constants
consensus = "modz"
data_dir = pathlib.Path("../1.generate-profiles/data/consensus")
output_dir = pathlib.Path("results")
output_dir.mkdir(exist_ok=True)
output_file = pathlib.Path(f"{output_dir}/cell_line_holdout_robustness_results_{consensus}.tsv.gz")

cell_lines = ["A549", "ES2", "HCC44"]
cell_health_metadata = ["Metadata_profile_id", "Metadata_pert_name", "Metadata_cell_line"]
shuffle_types = [True, False]


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
            max_iter=2000,
            tol=1e-3
        )
        
    )]
)

y_transform = "raw"
scoring = "r2"
decision_function = False
model_type = "Regression"


# In[5]:


# Load X data
x_file = pathlib.Path(f"{data_dir}/cell_painting_{consensus}.tsv.gz")
x_df = pd.read_csv(x_file, sep="\t")

print(x_df.shape)
x_df.head(2)


# In[6]:


# Load Y data
y_file = pathlib.Path(f"{data_dir}/cell_health_{consensus}.tsv.gz")
y_df = pd.read_csv(y_file, sep="\t")

# Setup target variable names and assert that there are 70 of them
cell_health_targets = y_df.drop(cell_health_metadata, axis="columns").columns.tolist()
assert len(cell_health_targets) == 70

print(y_df.shape)
y_df.head(2)


# In[7]:


regression_results_list = []
for cell_line in cell_lines:
    # Split into training and testing sets
    train_x_df = process_cell_line_train_test(x_df, cell_line, cell_health_metadata, test_set=False)
    train_y_df = process_cell_line_train_test(y_df, cell_line, cell_health_metadata, test_set=False)
    test_x_df = process_cell_line_train_test(x_df, cell_line, cell_health_metadata, test_set=True)
    test_y_df = process_cell_line_train_test(y_df, cell_line, cell_health_metadata, test_set=True)

    for cell_health_target in cell_health_targets:
        for shuffle_type in shuffle_types:
            print(f"Now Training Target: {cell_health_target}")
            print(f"[Class] Cell Line: {cell_line}; Shuffle: {shuffle_type}\n")

            # Initialize class
            chp = CellHealthPredict(
                x_df=train_x_df,
                y_df=train_y_df,
                parameters=parameters,
                estimator=estimator,
                n_folds=n_folds,
                cv_scoring=scoring,
                shuffle=shuffle_type
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
                cell_line=cell_line
            )

            # Testing performance metrics
            metric_mse_test, metric_rtwo_test, y_test_true, y_test_pred = chp.get_performance(
                x_test=test_x_df,
                y_test=test_y_df,
                decision_function=decision_function,
                return_y=True,
                data_fit_type="test",
                cell_line=cell_line
            )

            # Store results
            regression_results_list.append(pd.concat([metric_mse, metric_mse_test], axis='rows'))
            regression_results_list.append(pd.concat([metric_rtwo, metric_rtwo_test], axis='rows'))


# In[8]:


full_regression_results_df = pd.concat(regression_results_list).reset_index(drop=True)

print(full_regression_results_df.shape)
full_regression_results_df.head(3)


# In[9]:


# Save all results
full_regression_results_df.to_csv(output_file, sep="\t", index=False, compression="gzip")

