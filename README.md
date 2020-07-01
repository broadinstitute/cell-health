# Predicting cell health phenotypes using image-based morphology profiling

**Gregory P. Way, Maria Kost-Alimova, Tsukasa Shibue, William F. Harrington, Stanley Gill, Tim Becker, William C. Hahn, Anne E. Carpenter, Francisca Vazquez2, Shantanu Singh 2020**

## Summary

The following repository stores a complete analysis pipeline using Cell Painting data to predict readouts from several cell health assays.

![approach](https://raw.githubusercontent.com/broadinstitute/cell-health/master/media/approach.png)
This overview figure outlines the Cell Health assay, the Cell Painting assay, and our machine learning approach.

First, our collaborators developed a customized microscopy assay we call "Cell Health".
The Cell Health assay is comprised of two different reagent panels: Cell cycle and viability.
Together, these two panels use reagents which mark different cell health phenotypes.

| Assay/Dye | Phenotype | Panel |
| :-- | :-------- | :---------- |
| Caspase 3/7 | Apoptosis | Viability |
| DRAQ7 | Cell Death | Viability |
| CellROX | Reactive Oxygen Species | Viability
| EdU | Cellular Profileration | Cell Cycle |
| Hoechst | DNA Content | Cell Cycle |
| pH3 | Cell Division | Cell Cycle |
| gH2Ax | DNA Damage | Cell Cycle |

We hypothesized that we can use unbiased and high dimensional Cell Painting profiles to predict the readouts of each individual assay.

## Data

We collected Cell Painting measurements using CRISPR perturbations.
The experiment targeted 59 genes, which included 119 unique guides (~2 per gene), across 3 cell lines.
The cell lines included A549, ES2, and HCC44.

| Cell Line | Primary Site |
| :-------- | :----------- |
| A549 | Lung Cancer |
| ES2 | Ovarian Cancer |
| HCC44 | Lung Cancer |

![CRISPR Correlation](https://raw.githubusercontent.com/broadinstitute/cell-health/master/2.replicate-reproducibility/figures/guide_correlation.png)

About 60% of all CRISPR guides were reproducible.
This is consistent with previous genetic perturbations ([Rohban et al. 2017](https://doi.org/10.7554/eLife.24060)).
It is important to note that we are not actually interested in the CRISPR treatment specifically, but instead, just its corresponding readout in each cell health assay.

## Pipeline

The full analysis pipeline consists of the following steps:

| Order | Module | Description |
| :---- | :----- | :---------- |
| 0 | Download cell painting data | Retrieve single cell profiles archived on Figshare |
| 1 | Generate profiles | Generate and process cell painting and cell health assay readouts |
| 2 | Determine replicate reproducibility | Determine the extent to which the CRISPR perturbations result in reproducible signatures |
| 3 | Train machine learning models to predict cell health assays | Train and visualize regression models using cell painting data to predict cell health assay readouts |
| 4 | Apply the models | Apply the trained models to the Drug Repurposing Hub data to predict drug perturbation effect |
| 5 | Validate the models | Use orthogonal readouts to validate the Drug Repurposing Hub predictions |

Each analysis module should be run in order.
View each module for specific instructions on how to reproduce results.

[`analysis-pipeline.sh`](analysis-pipeline.sh) stores information on how to reproduce all analysis modules.

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

## Machine Learning Approach

We performed the following approach:

1. Split data into 85% training and 15% test sets.
2. Normalized data using the EMPTY controls in each plate (moderated z-score).
3. Selected optimal hyperparamters using 5-fold cross-validation.
4. Trained elastic net regression models to predict each of the 70 cell health assay readouts, independently.
5. Trained using shuffled data as well.
6. Report performance on training and test sets.

## Results

### Regression Model Performance

![Regression Model Performance](https://raw.githubusercontent.com/broadinstitute/cell-health/master/3.train/figures/regression/modz/regression_performance_figure_modz.png)

Initial results indicate that many of the cell health phenotypes can be predicted with high performance using our approach.
However, there are many cell line specific differences.

### Model Interpretation

![Model Interpretation](https://raw.githubusercontent.com/broadinstitute/cell-health/master/3.train/figures/coefficients/coefficient_summary_Real_modz_abs_mean.png)

Because we used a logistic regression classifier, we can readily interpret the output features.
These features were derived from CellProfiler and represent different measurements of cell morphology
Shown above is a summary of coefficients from all 70 cell health models.
We observed that each contribute to classifying various facets of cell health.
Many different categories of cell morphology features contribute to cell health predictions.

### Application to Drug Repurposing Hub

We applied the trained models to Cell Painting data from the Drug Repurposing Hub ([Corsello et al. 2017](https://doi.org/10.1038/nm.4306)).
These data represent ~1,500 compound perturbations in ~6 dose points in A549 cells.

We applied all predictions and present them in an easy-to-navigate webapp at https://broad.io/cell-health

![lincs](https://raw.githubusercontent.com/broadinstitute/cell-health/master/4.apply/figures/lincs_main_figure_4.png)

Collapsing the Drug Repurposing Hub Cell Painting data into UMAP coordinates, we observed many associated Cell Health predictions.
For example, predicted G1 Cell Count and predicted ROS had clear gradients in UMAP space.
However, there is not exactly a 1-1 relationship.
The control proteasome inhibitors (DMSO and Bortezomib) are known to induce ROS, while PLK inhibitors are known to induce cell death by blocking mitosis entry.
A single PLK inhibitor (HMN-214) showed a strong dose relationship with predicted G1 count.
