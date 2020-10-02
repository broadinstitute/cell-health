#!/usr/bin/env python
# coding: utf-8

# ## Performing a sample titration analysis
# 
# Our goal is to determine the effect of sample size on model performance.
# 
# ### Approach
# 
# We will use the existing training/test sets to evaluate overall performance on the test set by systematically reducing the samples included in the training set.
# 
# 1. Randomly shuffle the index (aka the order) of the training set samples
# 2. Train cell health model using all samples in the training set
# 3. Remove the top ordered sample from the training set and retrain
# 4. Repeat step 3 after removing up to 50 samples from training
# 5. Repeat step 1 (different random shuffling 100 times)
# 6. For each training iteration, collect the test set R2 performance
# 
# In total, this results in 70 (models) x 25 (different sample sets) x 3 (different iterations) = 5,640 different model initializations

# In[9]:


import sys
import pathlib
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
output_dir = pathlib.Path("results")
output_dir.mkdir(exist_ok=True)
output_file = pathlib.Path(f"{output_dir}/sample_titration_robustness_results_{consensus}.tsv.gz")
output_dropped_file = pathlib.Path(f"{output_dir}/sample_titration_samples_dropped_small_{consensus}.tsv.gz")

num_iterations = 3
num_sample_titration = 25


# In[3]:


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


# In[4]:


# Load data
x_train_df, x_test_df, y_train_df, y_test_df = (
    load_train_test(data_dir=data_dir, drop_metadata=True, consensus=consensus)
)

# Setup target variable names and assert that there are 70 of them
cell_health_targets = y_train_df.columns.tolist()
assert len(cell_health_targets) == 70


# In[5]:


regression_results_list = []
samples_dropped = []
for iteration in range(0, num_iterations):
    x_sample_df = x_train_df.sample(frac=1)
    for drop_sample_high in range(1, num_sample_titration+1):
        drop_samples = x_sample_df.iloc[range(0, drop_sample_high), :].index.tolist()
        
        x_train_subset_df = x_sample_df.drop(drop_samples, axis="index")
        y_train_subset_df = y_train_df.reindex(x_train_subset_df.index, axis="index")
        
        # Store knowledge of which samples were dropped per iteration
        drop_sampels_str = ";".join(drop_samples)
        samples_dropped.append(pd.Series([iteration, drop_sample_high, drop_sampels_str]))
        
        for cell_health_target in cell_health_targets:
            print(f"Now Training Target: {cell_health_target}")
            print(f"[Class] Titration: {len(drop_samples)}; Iteration: {iteration}\n")
            # Initialize class
            chp = CellHealthPredict(
                x_df=x_train_subset_df,
                y_df=y_train_subset_df,
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
                x_test=x_test_df,
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
                    num_samples_dropped=drop_sample_high,
                    iteration=iteration
                )
            )
            
            rtwo_df = (
                pd.concat([metric_rtwo, metric_rtwo_test], axis='rows')
                .reset_index(drop=True)
                .assign(
                    num_samples_dropped=drop_sample_high,
                    iteration=iteration
                )
            )
            
            regression_results_list.append(mse_df)
            regression_results_list.append(rtwo_df)


# In[6]:


# Compile the full regression results
full_regression_results_df = pd.concat(regression_results_list).reset_index(drop=True)

print(full_regression_results_df.shape)
full_regression_results_df.head(3)


# In[7]:


# Compile the samples dropped
samples_dropped_df = (
    pd.concat([pd.DataFrame(x).transpose() for x in samples_dropped])
    .reset_index(drop=True)
)
samples_dropped_df.columns = ["iteration", "num_dropped", "samples_dropped"]

print(samples_dropped_df.shape)
samples_dropped_df.head()


# In[8]:


# Save all results
full_regression_results_df.to_csv(output_file, sep="\t", index=False, compression="gzip")
samples_dropped_df.to_csv(output_dropped_file, sep="\t", index=False, compression="gzip")

