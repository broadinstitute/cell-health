# Download Single-Cell Cell Painting Profiles

**Gregory Way, 2019**

In this module, we download all single-cell Cell Painting profiles generated for the Cell Health project.
All of the single cell profiles are publicly available and hosted by [NIH Figshare](https://nih.figshare.com/).

The files are located at https://nih.figshare.com/account/articles/9995672.

## Recommendation

The size of the single cell profile data is ~130GB.

There is no need to download the single cell profiles.
We provide documentation here to access, if there are specific questions that only the single cell data can answer.

Instead, all processed data used for all downstream analyses is located in [1.generate-profiles/data/profiles/CRISPR_PILOT_B1/](https://github.com/broadinstitute/cell-health/tree/master/1.generate-profiles/data/profiles/CRISPR_PILOT_B1).

## Download

To download the single cell profiles, perform the following in your terminal:

```bash
# Activate environment
conda activate cell-health

# Navigate to appropriate directory
cd 0.download-data

# Run download script
python download.py
```

## Upload

**Note:** The upload script provided (`upload.py`) was not actually executed in this environment.
We previously executed the script in a private Amazon ec2 instance with access to the Carpenter-Lab S3 bucket.
This bucket stores CellProfiler output for the Cell Health project.

Nevertheless, to execute the upload call, we performed the following:

```bash
# First, generate a figshare API token from https://nih.figshare.com/account/applications

# Activate environment
conda activate cell-health

# Perform upload
python upload.py --token <INSERT TOKEN HERE>
```
