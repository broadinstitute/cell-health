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


# ## Setup Cross Validation

# In[4]:


alphas = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
l1_ratios = [0.1, 0.12, 0.14, 0.16, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9]
n_folds = 5


# In[5]:


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


# In[6]:


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


# In[7]:


# Y labels and transform instructions
cell_health_targets = y_train_df.columns.tolist()
y_transforms = ["raw", "zero-one", "binarize"]


# ## Train Models

# In[8]:


get_ipython().run_cell_magic('time', '', '\ncv_results_list = []\nroc_results_list = []\npr_results_list = []\nall_coefs_list = []\nall_y_label_list = []\nregression_results_list = []\ncannot_fit_list = []\n\nfor cell_health_target in cell_health_targets:\n    print(cell_health_target)\n    \n    for y_transform in y_transforms:\n\n        if y_transform == "binarize":\n            estimator = estimator_classifier\n            scoring = "roc_auc"\n            parameters = clf_parameters\n            decision_function = True\n        else:\n            estimator = estimator_regressor\n            scoring = "r2"\n            parameters = regression_parameters\n            decision_function = False\n    \n        for shuffle_type in [True, False]:\n\n            # Initialize class\n            chp = CellHealthPredict(\n                x_df=x_train_df,\n                y_df=y_train_df,\n                parameters=parameters,\n                estimator=estimator,\n                n_folds=n_folds,\n                cv_scoring=scoring,\n                shuffle=shuffle_type\n            )\n\n            # Fit model\n            is_fit = chp.fit_cell_health_target(cell_health_target,\n                                                y_transform=y_transform,\n                                                binarize_fit="median")\n            \n            if not is_fit:\n                cannot_fit_list.append([cell_health_target, y_transform, shuffle_type])\n                continue\n\n            # Get performance metrics for training and testing\n            metric_a, metric_b, y_true, y_pred = chp.get_performance(\n                decision_function=decision_function,\n                return_y=True,\n                binarize_fit="median"\n            )\n            metric_test_a, metric_test_b, y_test_true, y_test_pred = chp.get_performance(\n                x_test=x_test_df,\n                y_test=y_test_df,\n                decision_function=decision_function,\n                return_y=True,\n                binarize_fit="median"\n            )\n\n            # Combine training and testing results\n            if y_transform == "binarize":\n                roc_results_list.append(pd.concat([metric_a, metric_test_a], axis=\'rows\'))\n                pr_results_list.append(pd.concat([metric_b, metric_test_b], axis=\'rows\'))\n            else:\n                regression_results_list.append(pd.concat([metric_a, metric_test_a], axis=\'rows\'))\n                regression_results_list.append(pd.concat([metric_b, metric_test_b], axis=\'rows\'))\n\n            # Save cross validation results\n            cv_results_list.append(chp.get_cv_results())\n\n            # Save the model coefficients\n            model_file = "cell_health_target_{}_shuffle_{}_transform_{}.joblib".format(\n                cell_health_target, shuffle_type, y_transform\n            )\n            model_file = os.path.join("models", model_file)\n            coef_df = chp.get_coefficients(save_model=True, model_file=model_file)\n            all_coefs_list.append(coef_df)\n        \n            # Store y predictions recoded values\n            all_y_label_list.append(pd.concat([y_true, y_test_true, y_pred, y_test_pred]))')


# In[9]:


# Acquire output metrics
full_cv_df = pd.concat(cv_results_list).reset_index(drop=True)
full_regression_results_df = pd.concat(regression_results_list).reset_index(drop=True)
full_roc_df = pd.concat(roc_results_list).reset_index(drop=True)
full_pr_df = pd.concat(pr_results_list).reset_index(drop=True)
full_coef_df = pd.concat(all_coefs_list).reset_index(drop=True)
full_y_df = pd.concat(all_y_label_list).reset_index(drop=True)


# In[10]:


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

