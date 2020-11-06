#!/bin/bash
#
# Gregory Way 2020
# Predicting Cell Health
# Testing machine learning robustness
#
# 1) Assessing model performance on held out cell lines
# 2) Assessing model performance in a simulated situation where we have fewer samples
# 3) Assess model performance when feature groups are systematically heldout from analysis

set -e

# Step 0: Convert all notebooks to scripts
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts/nbconverted \
        *.ipynb

# Step 1: Cell line holdout analysis
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 0.cell-line-holdout.ipynb

# Step 2: Sample titration analysis
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000000 \
        --execute 1.sample-titration.ipynb

# Step 3: Feature group removal
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 2.feature-group-subsets.ipynb

# Step 4: Visualize results
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 3.visualize-robustness.ipynb
