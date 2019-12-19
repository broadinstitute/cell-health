#!/bin/bash

# Download all sqlite files from https://doi.org/10.35092/yhjc.9995672.v4
python download.py

# Confirm download integrity
md5sum -c md5sums.txt

