# Training Models to Predict Cell Health Phenotypes

**Gregory Way, 2020**

In the following module, we test the ability of Cell Painting to predict cell health readouts.

The `3.train` analysis module trains all models to predict cell health readouts.
The module stores data, software, results, figures, interpretations, and computational infrastructure to reproduce the analysis.

## Data

Both the Cell Painting and Cell Health assays were acquired on the same cell lines (A549, ES2, and HCC44) in the same perturbations.
The assays were acquired from a CRISPR experiment in which 59 genes were knocked out with, in most cases, 2 different CRISPR guides.

## Approach

First, we aggregated data at the _treatment_ level (aggregating replicate wells together) and matched cell health readouts (y) with Cell Painting profiles (X).
We randomly split the data (n = 357) into 85% training and 15% testing sets.

We then independently trained regression models to predict each of the 70 cell health readouts.
We also trained models again using randomly permuted Cell Painting data.

## Results

Many cell health readouts can be predicted with high accuracy in a testing set.

![Regression Model Performance](https://raw.githubusercontent.com/broadinstitute/cell-health/master/3.train/figures/regression/modz/regression_performance_figure_modz.png)

Some cell health readouts can be predicted exceptionally well, but others cannot be captured with our current approach.

![Individual Model Performance](https://raw.githubusercontent.com/broadinstitute/cell-health/master/3.train/figures/regression/modz/supplementary_figure_example_distributions.png)

Cell painting can be used to predict cell health readouts better than models trained with random.

![performance summary](https://raw.githubusercontent.com/broadinstitute/cell-health/master/3.train/figures/regression/modz/r_squared_comparison_scatter_modz.png)

More specific results exploring performance for all models and all cell health variables are presented in the `figures/` folder.
