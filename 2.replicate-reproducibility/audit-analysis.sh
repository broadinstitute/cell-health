#!/bin/bash
#
# Gregory Way 2019
# Predicting Cell Health
#
# Audit Analysis Pipeline

set -e

# Step 0 - Convert all notebooks to scripts
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted *.ipynb

# Step 1 - Analyze audits
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 0.analyze-audit.ipynb

# Step 2 - Visualize platemaps
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=ir \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.visualize-platemaps.ipynb

# Step 3 - Build morpheus gct files for heatmap generation
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 2.build-morpheus-gct.ipynb
