#!/bin/bash

##############################################
# Training Machine Learning Classifiers to Predict Cell Health Outcomes
#
# Gregory Way, 2019
##############################################

# Step 0: Convert all notebooks to scripts
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts/nbconverted \
        *.ipynb

# Step 1: Stratify data into training and testing sets
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 0.stratify-data.ipynb
