"""
Upload Single Cell Cell Painting Profiles
Gregory Way, 2019

We previously uploaded all single cell Cell Painting profiles to NIH Figshare.
See https://nih.figshare.com/account/articles/9995672 for more details.

The data are described in `0.download-data/download.py`.

**Important Note**

We previously executed this file, and running it will not work from a direct clone!
We executed this file on a private AWS ec2 instance with access to our CellProfiler S3 bucket.
"""

import argparse
from pycytominer.upload import figshare_upload

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--token", help="Figshare API Token")
args = parser.parse_args()

token = args.token
chunk_size = 1048576
article_id = "9995672"
backend_dir = "/home/ubuntu/bucket/projects/2015_07_01_Cell_Health_Vazquez_Cancer_Broad/workspace/backend/CRISPR_PILOT_B1"

for plate in os.listdir(backend_dir):
    plate_file = os.path.join(backend_dir, plate, "{}.sqlite".format(plate))
    figshare_upload(
        token=token,
        file_name=plate_file,
        append=True,
        article_id=article_id,
        chunk_size=chunk_size,
    )
