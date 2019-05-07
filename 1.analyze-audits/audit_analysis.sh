#!/bin/bash
#
# Gregory Way 2019
# Predicting Cell Health
#
# Audit Analysis Pipeline

# Step 0 - Convert all notebooks to scripts
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted *.ipynb

# Step 1 - Analyze audites
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute analyze-audit.ipynb
