#!/bin/bash
#
# Gregory Way 2020
# Predicting Cell Health
# Testing for correlations in A549 Viability Measurements
#
# I compare viability measurements between:
# 1) Cell Health Model Predictions (Number of Objects)
# 2) A549 PRISM Viability Screen

set -e

# Step 0: Convert all notebooks to scripts
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts/nbconverted \
        *.ipynb

# Step 1: Process Dependency Map Data
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 0.process-depmap-data.ipynb

# Step 2: Merge with cell health model and output results
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.validate.ipynb
