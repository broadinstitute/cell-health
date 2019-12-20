"""
Generate morphological profiles from cell painting data using pycytominer.

Do not perform feature selection since downstream analysis include feature selection
"""

import os
import numpy as np
import pandas as pd
import argparse
import multiprocessing
from joblib import Parallel, delayed

from pycytominer.aggregate import AggregateProfiles
from pycytominer import annotate, normalize, feature_select, audit


def get_profiles(plate, backend_dir, metadata_dir, barcode_platemap_df):
    """
    Apply all profiling steps for a given plate. To be applied in parallel.

    Output:
    Will write a series of processed files to disk
    """
    print("Processing {}.....".format(plate))
    plate_dir = os.path.join(backend_dir, plate)
    sqlite_file = "sqlite:////{}/{}.sqlite".format(plate_dir, plate)

    # Load specific platemap
    platemap = barcode_platemap_df.query(
        "Assay_Plate_Barcode == @plate"
    ).Plate_Map_Name.values[0]
    platemap_file = os.path.join(metadata_dir, "platemap", "{}.csv".format(platemap))
    platemap_df = pd.read_csv(platemap_file)

    # Prepare sql file for processing
    ap = AggregateProfiles(
        sqlite_file, strata=["Image_Metadata_Plate", "Image_Metadata_Well"]
    )

    # Count cells and output
    cell_count_file = os.path.join("results", "{}_cell_count.tsv".format(plate))
    cell_count_df = ap.count_cells()
    cell_count_df = cell_count_df.merge(
        platemap_df, left_on="Image_Metadata_Well", right_on="well_position"
    ).drop(["WellRow", "WellCol", "well_position"], axis="columns")
    cell_count_df.to_csv(cell_count_file, sep="\t", index=False)

    # Begin processing profiles
    output_dir = os.path.join("data", "profiles", batch, plate)
    os.makedirs(output_dir, exist_ok=True)

    # Aggregate single cells into well profiles
    out_file = os.path.join(output_dir, "{}.csv.gz".format(plate))
    ap.aggregate_profiles(output_file=out_file, compression="gzip")

    # Annotate Profiles
    anno_file = os.path.join(output_dir, "{}_augmented.csv.gz".format(plate))
    annotate(
        profiles=out_file,
        platemap=platemap_df,
        join_on=["Metadata_well_position", "Image_Metadata_Well"],
        output_file=anno_file,
        compression="gzip",
    )

    # Normalize Profiles
    norm_file = os.path.join(output_dir, "{}_normalized.csv.gz".format(plate))
    normalize(
        profiles=anno_file,
        features="infer",
        samples="Metadata_pert_name == 'EMPTY'",
        output_file=norm_file,
        compression="gzip",
    )

    # Perform feature selection (just drop columns with high number of missingness)
    feat_file = os.path.join(
        output_dir, "{}_normalized_feature_select.csv.gz".format(plate)
    )
    feature_select(
        profiles=norm_file,
        features="infer",
        samples="none",
        operation=["drop_na_columns", "blacklist"],
        output_file=feat_file,
        compression="gzip",
    )

    # Perform audits
    profile_df = pd.read_csv(feat_file).drop(
        ["Image_Metadata_Well", "Image_Metadata_Plate"], axis="columns"
    )

    # Audit guide replicability
    audit_file = os.path.join("results", "{}_audit_guide.csv".format(plate))
    audit(
        profiles=profile_df,
        audit_groups=["Metadata_pert_name", "Metadata_gene_name", "Metadata_cell_line"],
        iterations=10,
        output_file=audit_file,
    )

    # Audit gene replicability
    audit_file = os.path.join("results", "{}_audit_gene.csv".format(plate))
    audit(
        profiles=profile_df,
        audit_groups=["Metadata_gene_name", "Metadata_cell_line"],
        iterations=10,
        output_file=audit_file,
    )


parser = argparse.ArgumentParser()
parser.add_argument(
    "-p",
    "--parallel",
    action="store_true",
    help="decision to perform computation in parallel or not",
)
args = parser.parse_args()

parallel = args.parallel

num_cores = multiprocessing.cpu_count() - 1

batch = "CRISPR_PILOT_B1"
bucket_dir = os.path.join(
    "/home",
    "ubuntu",
    "bucket",
    "projects",
    "2015_07_01_Cell_Health_Vazquez_Cancer_Broad",
    "workspace",
)

backend_dir = os.path.join(bucket_dir, "backend", batch)
metadata_dir = os.path.join(bucket_dir, "metadata", batch)

# Load Barcode Platemap
barcode_platemap_file = os.path.join(metadata_dir, "barcode_platemap.csv")
barcode_platemap_df = pd.read_csv(barcode_platemap_file)

# Perform analysis for each plate
if __name__ == "__main__":
    all_plates = os.listdir(backend_dir)

    if parallel:
        Parallel(n_jobs=num_cores)(
            delayed(get_profiles)(
                plate=x,
                backend_dir=backend_dir,
                metadata_dir=metadata_dir,
                barcode_platemap_df=barcode_platemap_df,
            )
            for x in all_plates
        )

    else:
        for plate in all_plates:
            get_profiles(
                plate=plate,
                backend_dir=backend_dir,
                metadata_dir=metadata_dir,
                barcode_platemap_df=barcode_platemap_df,
            )
