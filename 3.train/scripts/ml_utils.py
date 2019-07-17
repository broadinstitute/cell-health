"""
Methods and functions to facilitate machine learning analyses of cell health data
"""

import os
import pandas as pd


def load_train_test(data_dir="data", drop_metadata=False, output_metadata_only=False):
    """
    Load training and testing data

    Arguments:
    data_dir - a string indicating the location of the files
    drop_metadata - boolean if the metadata should be striped before input
                    [default: False]
    output_metadata_only - boolean if metadata columns are only output
                           [default: False]

    Output:
    A list storing training and testing x and y matrices
    """

    x_train_file = os.path.join(data_dir, "x_train.tsv.gz")
    y_train_file = os.path.join(data_dir, "y_train.tsv.gz")
    x_test_file = os.path.join(data_dir, "x_test.tsv.gz")
    y_test_file = os.path.join(data_dir, "y_test.tsv.gz")

    x_train_df = pd.read_csv(x_train_file, sep="\t")
    y_train_df = pd.read_csv(y_train_file, sep="\t")
    x_test_df = pd.read_csv(x_test_file, sep="\t")
    y_test_df = pd.read_csv(y_test_file, sep="\t")

    x_train_metadata = x_train_df.columns.str.startswith("Metadata_")
    x_test_metadata = x_test_df.columns.str.startswith("Metadata_")

    y_metadata = ["guide", "cell_id"]

    if drop_metadata:
        x_train_df = x_train_df.loc[:, ~x_train_metadata]
        x_test_df = x_test_df.loc[:, ~x_test_metadata]

        y_train_df = y_train_df.drop(y_metadata, axis="columns")
        y_test_df = y_test_df.drop(y_metadata, axis="columns")

    elif output_metadata_only:
        x_train_df = x_train_df.loc[:, x_train_metadata]
        x_test_df = x_test_df.loc[:, x_test_metadata]

        y_train_df = y_train_df.loc[:, y_metadata]
        y_test_df = y_test_df.loc[:, y_metadata]

    return x_train_df, x_test_df, y_train_df, y_test_df
