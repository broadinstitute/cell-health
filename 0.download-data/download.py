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
import requests


def download_sqllite_file(filename, url):
    with requests.get(url, stream=True) as sql_request:
        sql_request.raise_for_status()
        with open(filename, 'wb') as sql_fh:
            for chunk in sql_request.iter_content(chunk_size=819200000):
                if chunk:
                    assert isinstance(chunk, object)
                    sql_fh.write(chunk)

file_info = {
    "SQ00014610": "18028784",
    "SQ00014611": "18508583",
    "SQ00014612": "18505937",
    "SQ00014613": "18506036",
    "SQ00014614": "18031619",
    "SQ00014615": "18506108",
    "SQ00014616": "18506912",
    "SQ00014617": "18508316",
    "SQ00014618": "18508421",
}

download_dir = "data"
os.makedirs(download_dir, exist_ok=True)

for plate in file_info:
    figshare_id = file_info[plate]
    filename = os.path.join(download_dir, "{}.sqlite".format(plate))
    if os.path.exists(filename):
        continue
    print("Now downloading... {}".format(filename))
    url = "https://nih.figshare.com/ndownloader/files/{}".format(figshare_id)
    download_sqllite_file(filename, url)
    print("Done...\n\n")
