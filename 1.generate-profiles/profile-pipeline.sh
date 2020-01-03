#!/bin/bash
#
# Gregory Way 2019
# Predicting Cell Health
# Pipeline to Process Data
#
# 1) Aggregate and Normalize Single Cell Data
# 2) Curate and Normalize Cell Health Assay Data

set -e

# Step 0: Convert all notebooks to scripts
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts/nbconverted \
        *.ipynb

# Step 1: Use pycytominer to generate and process all single cell profiles
python generate_profiles.py

# Step 2: Process cell health assay data (these are the labels for training models)
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 0.process-labels.ipynb

# Step 3: Normalize cell health assay data and generate .gct files for heatmaps
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.normalize-labels.ipynb

# Step 4: Build consensus signatures for downstream processing
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 2.build-consensus-signatures.ipynb
