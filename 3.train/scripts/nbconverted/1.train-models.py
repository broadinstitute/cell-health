#!/usr/bin/env python
# coding: utf-8

# # Train Classifier Models to Predict Cell Health Phenotypes
# 
# **Gregory Way, 2019**

# In[1]:


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.linear_model import SGDClassifier, ElasticNet
from sklearn.pipeline import Pipeline

from scripts.ml_utils import load_train_test, CellHealthPredict


# In[2]:


np.random.seed(123)


# ## Load Data

# In[3]:


x_train_df, x_test_df, y_train_df, y_test_df = load_train_test(drop_metadata=True)
x_meta_train_df, x_meta_test_df, y_meta_train_df, y_meta_test_df = load_train_test(output_metadata_only=True)


# In[4]:


cell_lines = list(set(x_meta_train_df.Metadata_cell_line))
cell_lines


# In[5]:


print(x_train_df.shape)
x_train_df.head(3)


# ## Setup Cross Validation

# In[6]:


alphas = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
l1_ratios = [0.1, 0.12, 0.14, 0.16, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9]
n_folds = 5


# In[7]:


regression_parameters = {
    'regress__alpha': alphas,
    'regress__l1_ratio': l1_ratios
}

clf_parameters = {
    'classify__loss': ['log'],
    'classify__penalty': ['elasticnet'],
    'classify__alpha': alphas,
    'classify__l1_ratio': l1_ratios
}


# In[8]:


estimator_regressor = Pipeline(
    steps=[(
        "regress",
        ElasticNet(
            random_state=42,
            max_iter=2000,
            tol=1e-3
        )
        
    )]
)

estimator_classifier = Pipeline(
    steps=[(
        'classify',
        SGDClassifier(
            random_state=42,
            class_weight='balanced',
            max_iter=2000,
            shuffle=True,
            tol=1e-3
        )
    )]
)


# In[9]:


# Y labels and transform instructions
cell_health_targets = y_train_df.columns.tolist()
y_transforms = ["raw", "binarize"]


# ## Train Models

# In[10]:


get_ipython().run_cell_magic('time', '', '\ncv_results_list = []\nroc_results_list = []\npr_results_list = []\nall_coefs_list = []\nall_y_label_list = []\nregression_results_list = []\ncannot_fit_list = []\n\nfor cell_health_target in cell_health_targets:\n    print(cell_health_target)\n    \n    for y_transform in y_transforms:\n\n        if y_transform == "binarize":\n            estimator = estimator_classifier\n            scoring = "roc_auc"\n            parameters = clf_parameters\n            decision_function = True\n        else:\n            estimator = estimator_regressor\n            scoring = "r2"\n            parameters = regression_parameters\n            decision_function = False\n    \n        for shuffle_type in [True, False]:\n\n            # Initialize class\n            chp = CellHealthPredict(\n                x_df=x_train_df,\n                y_df=y_train_df,\n                parameters=parameters,\n                estimator=estimator,\n                n_folds=n_folds,\n                cv_scoring=scoring,\n                shuffle=shuffle_type\n            )\n\n            # Fit model\n            is_fit = chp.fit_cell_health_target(cell_health_target,\n                                                y_transform=y_transform,\n                                                binarize_fit="median")\n            \n            if not is_fit:\n                cannot_fit_list.append([cell_health_target, y_transform, shuffle_type])\n                continue\n\n            # Get performance metrics for training and testing\n            metric_a, metric_b, y_true, y_pred = chp.get_performance(\n                decision_function=decision_function,\n                return_y=True,\n                binarize_fit="median"\n            )\n            metric_test_a, metric_test_b, y_test_true, y_test_pred = chp.get_performance(\n                x_test=x_test_df,\n                y_test=y_test_df,\n                decision_function=decision_function,\n                return_y=True,\n                binarize_fit="median",\n                data_fit_type="test",\n            )\n\n            # Get Cell Line Specific Performance\n            cell_line_metrics_a = []\n            cell_line_metrics_b = []\n            for cell_line in cell_lines:\n                meta_train_subset_df = x_meta_train_df.query("Metadata_cell_line == @cell_line")\n                meta_test_subset_df = x_meta_test_df.query("Metadata_cell_line == @cell_line")\n\n                # Get Cell Line Specific Training Performance\n                x_cell_line_df = x_train_df.reindex(meta_train_subset_df.index, axis="rows")\n                y_cell_line_df = y_train_df.reindex(meta_train_subset_df.index, axis="rows")\n\n                metric_cell_train_a, metric_cell_train_b, y_cell_train_true, y_cell_train_pred = (\n                    chp.get_performance(\n                        x_test=x_cell_line_df,\n                        y_test=y_cell_line_df,\n                        decision_function=decision_function,\n                        return_y=True,\n                        binarize_fit="median",\n                        cell_line=cell_line\n                    )\n                )\n\n                # Get Cell Line Specific Test Performance\n                x_cell_line_df = x_test_df.reindex(meta_test_subset_df.index, axis="rows")\n                y_cell_line_df = y_test_df.reindex(meta_test_subset_df.index, axis="rows")\n\n                metric_cell_test_a, metric_cell_test_b, y_cell_test_true, y_cell_test_pred = (\n                    chp.get_performance(\n                        x_test=x_cell_line_df,\n                        y_test=y_cell_line_df,\n                        decision_function=decision_function,\n                        return_y=True,\n                        binarize_fit="median",\n                        data_fit_type="test",\n                        cell_line=cell_line\n                    )\n                )\n\n                cell_line_metrics_a += [metric_cell_train_a, metric_cell_test_a]\n                cell_line_metrics_b += [metric_cell_train_b, metric_cell_test_b]\n\n            # Combine training and testing results\n            if y_transform == "binarize":\n                roc_results_list.append(pd.concat([metric_a, metric_test_a], axis=\'rows\'))\n                roc_results_list.append(pd.concat(cell_line_metrics_a, axis="rows"))\n                pr_results_list.append(pd.concat([metric_b, metric_test_b], axis=\'rows\'))\n                pr_results_list.append(pd.concat(cell_line_metrics_b, axis="rows"))\n            else:\n                regression_results_list.append(pd.concat([metric_a, metric_test_a], axis=\'rows\'))\n                regression_results_list.append(pd.concat([metric_b, metric_test_b], axis=\'rows\'))\n                regression_results_list.append(pd.concat(cell_line_metrics_a, axis=\'rows\'))\n                regression_results_list.append(pd.concat(cell_line_metrics_b, axis="rows"))\n\n            # Save cross validation results\n            cv_results_list.append(chp.get_cv_results())\n\n            # Save the model coefficients\n            model_file = "cell_health_target_{}_shuffle_{}_transform_{}.joblib".format(\n                cell_health_target, shuffle_type, y_transform\n            )\n            model_file = os.path.join("models", model_file)\n            coef_df = chp.get_coefficients(save_model=True, model_file=model_file)\n            all_coefs_list.append(coef_df)\n        \n            # Store y predictions recoded values\n            all_y_label_list.append(pd.concat([y_true, y_test_true, y_pred, y_test_pred]))')


# In[11]:


# Acquire output metrics
full_cv_df = pd.concat(cv_results_list).reset_index(drop=True)
full_regression_results_df = pd.concat(regression_results_list).reset_index(drop=True)
full_roc_df = pd.concat(roc_results_list).reset_index(drop=True)
full_pr_df = pd.concat(pr_results_list).reset_index(drop=True)
full_coef_df = pd.concat(all_coefs_list).reset_index(drop=True)
full_y_df = pd.concat(all_y_label_list).reset_index(drop=True)


# In[12]:


# Save all results
results_dir = "results"
os.makedirs(results_dir, exist_ok=True)

file = os.path.join(results_dir, "full_cell_health_cv_results.tsv.gz")
full_cv_df.to_csv(file, sep='\t', index=False)

file = os.path.join(results_dir, "full_cell_health_regression_results.tsv.gz")
full_regression_results_df.to_csv(file, sep='\t', index=False)

file = os.path.join(results_dir, "full_cell_health_roc_results.tsv.gz")
full_roc_df.to_csv(file, sep='\t', index=False)

file = os.path.join(results_dir, "full_cell_health_pr_results.tsv.gz")
full_pr_df.to_csv(file, sep='\t', index=False)

file = os.path.join(results_dir, "full_cell_health_coefficients.tsv.gz")
full_coef_df.to_csv(file, sep='\t', index=False)

file = os.path.join(results_dir, "full_cell_health_y_labels.tsv.gz")
full_y_df.to_csv(file, sep='\t', index=False)

