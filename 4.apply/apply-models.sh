#!/bin/bash

##############################################
# Applying trained models to Cell Painting Data from
# The Drug Repurposing Hub
#
# Gregory Way, 2020
##############################################

# Step 0: Convert all notebooks to scripts
jupyter nbconvert --to=script \
        --FilesWriter.build_directory=scripts/nbconverted \
        *.ipynb

# Step 1: Apply the models to the Drug Repurposing Hub
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 0.apply-models-repurposing.ipynb

# Step 2: Visualize results
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.visualize-repurposing

# Step 3: Fit dose curves
Rscript --vanilla 2.fit-dose.R

# Step 4: Visualize dose response curves
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 3.visualize-dose-response.ipynb

# Step 5: Generate summary figures
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 4.lincs-figures.ipynb
