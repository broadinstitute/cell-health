#!/bin/bash

##############################################
# Training Machine Learning Classifiers to Predict Cell Health Outcomes
#
# Gregory Way, 2020
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

# Step 2: Train all the machine learning models (classification and regression)
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.train-models.ipynb

# Step 3: Visualize distribution of features
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 2.visualize-binary-distribution.ipynb

# Step 4A: Visualize Performance (Classification)
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 3.visualize-classification-performance.ipynb

# Step 4B: Visualize Performance (Regression)
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 4.visualize-regression-performance.ipynb

# Step 5: Summarize model performance
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 5.visualize-performance-summary.ipynb

# Step 6: Apply models to full dataset
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 6.apply-models.ipynb

# Step 7: Visualize cell line specific performance
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 7.visualize-cell-line-performance.ipynb

# Step 8: Visualize model coefficients
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 8.visualize-coefficients.ipynb

# Step 9: Visualize example cells
# Note: This step requires the raw image data!
# jupyter nbconvert --to=html \
#         --FilesWriter.build_directory=scripts/html \
#         --ExecutePreprocessor.kernel_name=ir \
#         --ExecutePreprocessor.timeout=10000000 \
#         --execute 8.visualize-coefficients.ipynb
