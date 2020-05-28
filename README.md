# Predicting Cell Health with Morphological Profiles

**Gregory Way, Shantanu Singh, Anne Carpenter 2019**

## Summary

The following repository stores a complete analysis pipeline of using Cell Painting data to predict the results of several cell health assays.

Our collaborators previously collected a series of assays measuring different aspects of cell health.
The assays include staining with specific dyes to measure specific phenotypes.

| Assay/Dye | Phenotype |
| :-- | :-------- |
| Caspase 3/7 | Apoptosis |
| DRAQ7 | Cell Death |
| EdU | Cellular Profileration |
| Hoechst | DNA Content |
| Ph3 | Cell Division |
| gH2Ax | DNA Damage |

We hypothesized that we can use unbiased and high dimensional Cell Painting profiles to predict the readouts of each individual assay.

## Data

We collected Cell Painting measurements on a CRISPR experiment.
The experiment targeted 59 genes, which included 119 unique guides (~2 per gene), across 3 cell lines.
The cell lines included A549, ES2, and HCC44.

![CRISPR Correlation](https://raw.githubusercontent.com/broadinstitute/cell-health/master/2.replicate-reproducibility/figures/guide_correlation.png)

About 40% of all CRISPR guides were reproducible.
This is ok since we are not actually interested in the CRISPR treatment specifically, but instead, just its corresponding readout in each cell health assay.

## Pipeline

The full analysis pipeline consists of the following steps:

| Order | Module | Description |
| :---- | :----- | :---------- |
| 0 | Download cell painting data | Retrieve single cell profiles archived on Figshare |
| 1 | Generate profiles | Generate and process cell painting and cell health assay readouts |
| 2 | Determine replicate reproducibility | Determine the extent to which the CRISPR perturbations result in reproducible signatures |
| 3 | Train machine learning models to predict cell health assays | Train regression and classification models using cell painting data to predict cell health assay readouts |
| 4 | Apply the models | Apply the trained models to the Drug Repurposing Hub data to predict drug perturbation effect |

Each analysis module should be run in order.
View each module for specific instructions on how to reproduce results.

[`analysis-pipeline.sh`](analysis-pipeline.sh) stores information on how to reproduce all analysis modules.

## Machine Learning Approach

We performed the following approach:

1. Split data into 85% training and 15% test sets.
2. Normalized data by plate (z-score).
3. Selected optimal hyperparamters using 5-fold cross-validation
4. Trained elastic net regression models to predict each of the 70 cell health assay readouts, independently.
5. Trained using shuffled data as well.
6. Report performance on training and test sets.

We also trained logistic regression classifiers using the same approach above.

## Results

![Regression Model Performance](https://raw.githubusercontent.com/broadinstitute/cell-health/master/3.train/figures/regression/modz/regression_performance_figure_modz.png)

Initial results indicate that many of the cell health phenotypes can be predicted with our approach.
However, there are many cell line specific differences.

## Computational Environment

We use conda as a package manager.
To install conda see [instructions](https://docs.conda.io/en/latest/miniconda.html).

We recommend installing conda by downloading and executing the `.sh` file and accepting defaults.

To create the computational environment, run the following:

```sh
# Make sure the repo is cloned
conda env create --force --file environment.yml
conda activate cell-health
```
