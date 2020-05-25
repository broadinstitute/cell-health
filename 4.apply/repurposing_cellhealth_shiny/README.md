# Exploring Cell Health Predictions for Repurposing Hub Compounds

We present a shiny app to enable easy exploration of cell health predictions for 1,571 Drug Repurposing compound perturbations.
The app can be viewed at http://broad.io/cell-health-app.

Previously, we processed **a subset** of the [Broad Drug Repurposing Hub](https://clue.io/repurposing#home) compounds.
The processing scripts and data can be found [here](https://github.com/broadinstitute/lincs-cell-painting).
We use these processed data in this application.

The data are Cell Painting morphology profiles that are aggregated across single cells and biological replicates to form consensus perturbation signatures.
We test all perturbations across ~6 doses per compound.

We apply all [70 cell health machine learning models](https://github.com/broadinstitute/cell-health/blob/master/1.generate-profiles/data/labels/feature_mapping_annotated.csv) to the Drug Repurposing collection and present the results in the shiny app.
The following document provides instructions on how to use the app.

## Predicting Cell Health

This application represents one piece of a larger project that aims to determine how well certain cell health measurements can be predicted from a Cell Painting assay.
For more details about this larger project, see https://github.com/broadinstitute/cell-health.

## Getting Started

The app can be viewed at http://broad.io/cell-health-app.
There are two separate tabs, which are each useful for different explorations.
The `Model Explorer` tab enables quick exploration of cell health predictions across the profiled Drug Repurposing compounds.
The `Compound Explorer` tab focuses on specific user-specified compounds and multiple cell health model predictions.
We will describe specific options for each tab below:

## Model Explorer

This view focuses on three main plots.

1. Scatter Plot
    * The scatter plot will show either continuous model scores applied to the Drug Repurposing Hub data or a UMAP Of Drug Repurposing Hub consensus morphology profiles.
2. Model Ranking
    * All models were trained using external data (see https://github.com/broadinstitute/cell-health).
    * The bar chart shows performance of the cell health models in a holdout test set of the external data.
    * The performance of the models represents how we expect the models to perform on new data.
    * The bigger the bar, the more confidence we have in the subsequent predictions.
    * The color of the bar represents the dyes or antibodies used to generate the specific Cell Health measurement.
3. Dose Information
    * If a compound and cell health model score has dose support, then we are more confident in the models' predictions.
    * We show two dose plots with control points for the two axes given with the `Cell Health` scatter plot type.

The left side panel provides several options to interact with the compound perturbations and cell health predictions.
[Click here](https://github.com/broadinstitute/cell-health/blob/master/1.generate-profiles/data/labels/feature_mapping_annotated.csv) for a detailed description of each cell health measurement.

### Select Scatter Plot Type

There are two different scatter plot views available:

1. Cell Health [default]
    * Use this view to compare two separate cell health models and to isolate specific compounds of interest
    * This view represents continuous model prediction scores made for all cell health models
    * The color of the points represents the dose (ranked from low (1) to high (7))
2. UMAP
    * Use this view to determine how specific cell health model predictions and compounds are distributed according to global data structure
    * This view represents the consensus profiles transformed using UMAP ([McInnes et al. 2018](https://arxiv.org/abs/1802.03426))
    * The color of the points can be toggled, and represent continuous activation scores

In both cases, we highlight three control perturbations: DMSO (negative control), Bortezomib, and MG-132.
Bortezomib and MG-132 are both proteasome inhibitors and their impact on cell morphology is well characterized.
Note that dose information has been recoded to ["dose levels"](https://github.com/broadinstitute/lincs-cell-painting/tree/master/consensus#recoding-dose-information).

The compound selected in "Select a Compound" (see below) is highlighted in large triangles.
The controls are also highlighted by transparent grey shapes.

#### Scatter Plot Interaction

The user can click and draw a box around points for both scatter plot types.
Information about the selected points is displayed at the bottom of the page.
The information includes:

* Perturbation (pert_iname)
* Mechanism of Action (MOA)
* Genetic Target
* Clinical Phase
* Sample ID (Broad ID)
* Dose (Recoded)
* X and Y Axis Coordinates

For more details about the compounds see https://clue.io/repurposing-app

### Select Cell Health Variable for Y axis (UMAP Toggle)

The effect of this option depends on what is selected for the scatter plot type.
The option controls the Y axis with a `Cell Health` scatter plot type and the color option with a `UMAP` scatter plot type.
A detailed description of each cell health measurement can be found [here](https://github.com/broadinstitute/cell-health/blob/master/1.generate-profiles/data/labels/feature_mapping_annotated.csv).

### Select Cell Health Variable for X axis

Control the X axis with a `Cell Health` scatter plot type.
This option has no effect with a `UMAP` scatter plot type.

### Select a Compound

This option enables a user to select a specific compound to track in the scatter plot and dose views.
The form can accept copy and paste and will perform autocomplete.

### Remove Controls from "Click and Drag"

Check this box to remove DMSO, Bortezomib, and MG-132 from the table view.

## Compound Explorer

This tab focuses exclusively on a single compound at a time.
Select the specific combination of cell health measurements to visualize.
Compare your favorite compound of interest to controls and every other perturbation.

We also highlight which cell health models are being viewed with their associated expected performance.

### Select a Compound

Select which compound to focus on.
The form can accept copy and paste and will perform autocomplete.

### Check Models to Visualize

Use check boxes to select which cell health variables should be highlighted in the provided figures.
