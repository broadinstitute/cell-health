# Predicting cell health phenotypes using image-based morphology profiling

**2020**

**Gregory P. Way, Maria Kost-Alimova, Tsukasa Shibue, William F. Harrington, Stanley Gill, Tim Becker, William C. Hahn, Anne E. Carpenter, Francisca Vazquez, Shantanu Singh**

## Summary

The following repository stores a complete analysis pipeline using Cell Painting data to predict readouts from several cell health assays.

We first developed a customized microscopy assay we call "Cell Health".
The Cell Health assay is comprised of two different reagent panels: "Cell cycle" and "viability".
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

## Approach

This overview figure outlines the Cell Health assay, the Cell Painting assay, and our machine learning approach.

![approach](https://raw.githubusercontent.com/broadinstitute/cell-health/master/media/approach.png)

> Data processing and modeling approach.
> (a) 70 different readouts from the Cell Health imaging assay.
> The assay consists of two customized reagent panels, which use measurements from seven different targeted reagents and one channel based on digital phase contrast (DPC) imaging; shown are five toy examples to demonstrate that individual cells are isolated into subpopulations by various gating strategies to define the cell health readouts.
> (b) Example Cell Painting image across five channels, plus a merged representation across channels.
> The example image is cropped from a larger image and shows ES2 cells.
> Scale bars are 20 µm.
> Below are the steps applied in an image-based profiling pipeline, after image-based features have been extracted from each cell’s image.
> (c) Modeling approach where we fit 70 different regression models using CellProfiler features derived from Cell Painting images to predict Cell Health readouts.

## Data

### Access

All data are publicly available.

Single cell profiles are available at https://doi.org/10.35092/yhjc.9995672.
All processed data are available in this repository (see [`1.generate-profiles`](1.generate-profiles/)).

### Summary

We collected Cell Painting measurements using CRISPR perturbations.
The experiment targeted 59 genes, which included 119 unique guides (~2 per gene), across 3 cell lines.
The cell lines included A549, ES2, and HCC44.

| Cell Line | Primary Site |
| :-------- | :----------- |
| A549 | Lung Cancer |
| ES2 | Ovarian Cancer |
| HCC44 | Lung Cancer |

About 60% of all CRISPR guides were reproducible.
This is consistent with previous genetic perturbations ([Rohban et al. 2017](https://doi.org/10.7554/eLife.24060)).
It is important to note that we are not actually interested in the CRISPR treatment specifically, but instead, just its corresponding readout in each cell health assay.

![CRISPR Correlation](https://raw.githubusercontent.com/broadinstitute/cell-health/master/2.replicate-reproducibility/figures/guide_correlation.png)

> Median pairwise Pearson correlation of CRISPR guide replicate profiles (y axis) compared against Median pairwise Pearson correlation of CRISPR guides targeting the same gene or construct.
> We removed biological replicates when calculating the same-gene correlations.
> The three different cell lines (A549, ES2, and HCC44) are shown in different colors and in different facets of the figure.
> We generated the profiles by median aggregating CellProfiler measurements for all single cells within each well of a Cell Painting experiment (see Methods for more processing details).
> The text labels represent the proportion of gene and guide profiles with “strong phenotypes”.
> In other words, these profiles had replicate correlations greater than 95% of non-replicate pairwise Pearson correlations in the particular cell line.
> The dotted red line represents this 95% cutoff in the null distribution and the blue dotted line is y = x, which shows a strong consistency across CRISPR guide constructs.

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

### Computational Environment

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

Initial results indicate that many of the cell health phenotypes can be predicted with high performance using our approach.
However, there are many cell line specific differences.

![Regression Model Performance](https://raw.githubusercontent.com/broadinstitute/cell-health/master/3.train/figures/regression/modz/concise_regression_summary_modz.png)

> Test set model performance of predicting 70 cell health readouts with independent regression models.
> Performance for each phenotype is shown, sorted by decreasing R2 performance.
> The bars are colored based on the primary measurement metadata (see Supplementary Table S3), and they represent performance aggregated across the three cell lines.
> The points represent cell line specific performance. Points falling below -1 are truncated to -1 on the x axis. See 3.train for a full depiction.

### Model Interpretation

Because we used a logistic regression classifier, we can readily interpret the output features.
These features were derived from CellProfiler and represent different measurements of cell morphology
Shown above is a summary of coefficients from all 70 cell health models.
We observed that each contribute to classifying various facets of cell health.
Many different categories of cell morphology features contribute to cell health predictions.

![Model Interpretation](https://raw.githubusercontent.com/broadinstitute/cell-health/master/3.train/figures/coefficients/coefficient_summary_Real_modz_abs_mean.png)

> The importance of each class of Cell Painting features in predicting 70 cell health readouts.
> Each square represents the mean absolute value of model coefficients weighted by test set R2 across every model.
> The features are broken down by compartment (Cells, Cytoplasm, and Nuclei), channel (AGP, Nucleus, ER, Mito, Nucleolus/Cyto RNA), and feature group (AreaShape, Neighbors, Channel Colocalization, Texture, Radial Distribution, Intensity, and Granularity).
> For a complete description of all features, see the handbook: http://cellprofiler-manual.s3.amazonaws.com/CellProfiler-3.0.0/index.html
> Dark gray squares indicate “not applicable”, meaning either that there are no features in the class or the features did not survive an initial preprocessing step.
> Note that for improved visualization we multiplied the actual model coefficient value by 100.

### Application to Drug Repurposing Hub

We applied the trained models to Cell Painting data from the Drug Repurposing Hub ([Corsello et al. 2017](https://doi.org/10.1038/nm.4306)).
These data represent ~1,500 compound perturbations in ~6 dose points in A549 cells.

Collapsing the Drug Repurposing Hub Cell Painting data into UMAP coordinates, we observed many associated Cell Health predictions.
For example, predicted G1 Cell Count and predicted ROS had clear gradients in UMAP space.
However, there is not exactly a 1-1 relationship.
The control proteasome inhibitors (DMSO and Bortezomib) are known to induce ROS, while PLK inhibitors are known to induce cell death by blocking mitosis entry.
A single PLK inhibitor (HMN-214) showed a strong dose relationship with predicted G1 count.

![lincs](https://raw.githubusercontent.com/broadinstitute/cell-health/master/4.apply/figures/lincs_main_figure_4.png)

> Applying cell health models to Cell Painting data from The Drug Repurposing Hub.
> (a) We apply a Uniform Manifold Approximation (UMAP) to Drug Repurposing Hub consensus profiles of 1,571 compounds across 6 doses.
> The models were not trained using the Drug Repurposing Hub data.
> The point color represents the output of the cell health model trained to predict the number of cells in G1 phase (G1 cell count).
> (b) The same UMAP dimensions, but colored by the output of the Cell Health model trained to predict reactive oxygen species (ROS).
> (c) Drug Repurposing Hub profiles stratified by G1 cell count and ROS predictions.
> Bortezomib and MG-132 are proteasome inhibitors and are used as positive controls; DMSO is a negative control.
> We also highlight all PLK inhibitors in the dataset.
> (d) HMN-214 is an example of a PLK inhibitor that shows strong dose response for G1 cell count predictions.

We applied all predictions and present them in an easy-to-navigate webapp at https://broad.io/cell-health
