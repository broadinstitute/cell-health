"""
Download Single Cell Cell Painting Profiles
Gregory Way, 2019

These data were output from a CRISPR and Cell Painting experiment.
A microscope took pictures of all cells under several CRISPR perturbations.
The data were processed by a custom CellProfiler pipeline.
In the pipeline, we measured several morphology features for every single cell.
We then aggregated the output of CellProfiler using cytominer-database.
Cytominer-database compiles all single cell profiles into a single `sqlite` database.
We uploaed the `sqlite` files to NIH Figshare (see `0.download-data/upload.py`).
These data are publicly available

See https://nih.figshare.com/account/articles/9995672 for more details.

*Important Note*

These files are large (~130 GB Total).
There is no need to download and reprocess the data.
Processed data for all downstream analysis is available here:
    `1.generate-profiles/data/profiles/CRISPR_PILOT_B1/`.
"""

import os
from urllib.request import urlretrieve

file_info = {
    "SQ00014610": "18028784",
    "SQ00014611": "",
    "SQ00014612": "",
    "SQ00014613": "",
    "SQ00014614": "18031619",
    "SQ00014615": "",
    "SQ00014616": "",
    "SQ00014617": "",
    "SQ00014618": "",
}

download_dir = "data"
os.makedirs(download_dir, exist_ok=True)

for plate in file_info:
    figshare_id = file_info[plate]
    filename = os.path.join(download_dir, ".sqlite".format(plate))
    url = "https://nih.figshare.com/ndownloader/files/{}".format(figshare_id)
    urlretrieve(url, filename)
