"""
Generate morphological profiles from cell painting data using pycytominer.

Do not perform redundancy feature selection since downstream analysis uses models that induce sparsity.
"""

import os
import pandas as pd
import argparse
import multiprocessing

from pycytominer.aggregate import AggregateProfiles
from pycytominer import annotate, normalize, feature_select, audit


def get_profiles(plate, backend_dir, metadata_dir, barcode_platemap_df):
    """
    Apply all profiling steps for a given plate.

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
    output_dir = os.path.join("data", "profiles", plate)
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
        method="robustize",
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
        operation=["drop_na_columns", "blacklist", "variance_threshold"],
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


backend_dir = os.path.join("..", "0.download-data", "data")
metadata_dir = os.path.join("data", "metadata")

# Load Barcode Platemap
barcode_platemap_file = os.path.join(metadata_dir, "barcode_platemap.csv")
barcode_platemap_df = pd.read_csv(barcode_platemap_file)

# Perform analysis for each plate
if __name__ == "__main__":
    all_plates = [x.strip(".sqlite") for x in os.listdir(backend_dir) if x != ".DS_Store"]

    for plate in all_plates:
        get_profiles(
            plate=plate,
            backend_dir=backend_dir,
            metadata_dir=metadata_dir,
            barcode_platemap_df=barcode_platemap_df,
        )
