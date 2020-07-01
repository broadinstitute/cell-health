# Processing Data

**Gregory Way, 2020**

In this module, we present our pipeline for generating image-based profiles from Cell Painting data.
We also process Cell Health readouts.

## Generating image-based profiles

We primarily use the [pycytominer](https://github.com/cytomining/pycytominer) tool for data processing.

See `generate_profiles.py` for a complete description of our data processing pipeline.

Briefly, our pipeline is as follows:

| Step | Notes |
| :--- | :---- |
| Aggregate single cells | Operation: median |
| Annotate profiles | Merge platemaps with metadata |
| Normalize profiles | Operation: mad_robustize; using only EMPTY control wells |
| Feature select profiles | Operations: drop_na_columns, blacklist, variance_threshold, drop_outliers |
| Audit profiles | Determine quality of the data by pairwise replicate correlations |

## Processing cell health readouts

We also normalize the output cell health readouts from the Cell Health assay.
We simply take the z-score across features.

## Generating consensus signatures

We acquire consensus signatures for both Cell Painting and Cell Health assay readouts.
We generate two different types of consensus signatures: moderated z score (MODZ) and median consensus.

We use the MODZ operation in all downstream applications and interpretations.
MODZ was first introduced in [Subramanian et al., 2017](https://doi.org/10.1016/j.cell.2017.10.049) and we use the [pycytominer implementation](https://github.com/cytomining/pycytominer/blob/master/pycytominer/consensus.py).

This procedure results in a total of 357 profiles with matched Cell Painting and Cell Health data.

## Execution

To reprocess the profiles, execute the following command:

```bash
# Activate environment
conda activate cell-health

# Perform full profiling pipeline
# Note that step 6 of this pipeline is not currently executed,
# since raw images are required and not included in this repo.
python profile-pipeline.sh
```
