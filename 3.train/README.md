# Training Models to Predict Cell Health Phenotypes

**Gregory Way, 2019**

In the following module, we test the ability of Cell Painting ([Bray et al. 2016](https://doi.org/10.1038/nprot.2016.105)) to predict cell health assays.

Cell health assays represent 70 image based assays that use various dyes to capture important phenotypes in cell lines.
These phenotypes include proportion of cells in specific cell cycle stages, the kind of apoptosis initiated, cellular proliferation rates, reactive oxygen species (ROS) present, and many more.

The `3.train` analysis module trains all models to predict cell health readouts.
The module stores data, software, results, figures, interpretations, and computational infrastructure to reproduce the analysis.

## Data

Both the Cell Painting and Cell Health assays were acquired on the same cell lines (A549, ES2, and HCC44) in the same perturbations.
The assays were acquired from a CRISPR experiment in which 65 genes and controls were knocked out with, in most cases, 2 different CRISPR guides.

| Cell Line | Primary Site |
| :-------- | :----------- |
| A549 | Lung Cancer |
| ES2 | Ovarian Cancer |
| HCC44 | Lung Cancer |

## Approach

First, we aggregated data at the _treatment_ level (aggregating replicate wells together) and matched cell health readouts (y) with Cell Painting profiles (X).
We randomly split the data (n = 357) into 85% training and 15% testing sets.

We then independently trained three different models to predict each of the 70 cell health readouts.

| Model | Data Transformation |
| :-----| :------------------ |
| Classification | Binarize Cell Health Data by Median (after z-scoring) |
| Regression | Raw Cell Health Data (after z-scoring) |
| Regression | Zero-One Transform Cell Health Data (after z-scoring) |

We also trained all of these models again using randomly permuted cell painting data.

## Results

The cell health readouts can be predicted well by both classification and regression algorithms.

![performance summary](https://raw.githubusercontent.com/broadinstitute/cell-health/master/3.train/figures/performance_summary.png)

> Comparing regression mean squared error (y axis) to classification area under the receiver operating characteristic curve (AUROC) (x axis) in the test set.
The points represent the various cell health readouts colored by measurement category.
The shape of the point indicates if a collaborator (Maria Alimova) thought the measurements were particularly interesting.

Cell painting can be used to predict cell health readouts in regression models better than random.

![performance summary](https://raw.githubusercontent.com/broadinstitute/cell-health/master/3.train/figures/mse_comparison_scatter.png)

> Comparing permuted data mean squared error (y axis) to real data mean squared error (x axis) in raw cell health predictions.
Points above the red line indicate models that predicted cell health phenotypes better than models trained using permuted data.

More specific results exploring performance for all models and all cell health variables are presented in the `figures/` folder.
