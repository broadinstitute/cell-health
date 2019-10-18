#!/bin/bash
#
# Gregory Way 2019
# Predicting Cell Health
# Complete Analysis Pipeline

# Note that single cell profiles were uploaded to figshare and are available here:
# https://nih.figshare.com/account/articles/9995672
# To access this data, see the `0.download-data` module

# Step 1 - Generate Cell Painting profiles and normalize cell health assay readouts
cd 1.generate-profiles
bash profile-pipeline.sh

# Step 2 - Analyze replicate reproducibility
cd ..
cd 2.replicate-reproducibility
bash audit-analysis.sh

# Step 3 - Train machine learning models to predict cell health readouts
cd ..
cd 3.train
bash train-pipeline.sh

# Step 4 - Apply trained models to drug repurposing hub dataset
cd ..
cd 4.apply
bash apply-models.sh
