"""
Methods to visualize cell images.
The logic was adapted from:
https://github.com/jr0th/segmentation/blob/master/visualization/CellArt.ipynb
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import skimage.io
import skimage.exposure

from pycytominer.aggregate import AggregateProfiles


class grabPicture:
    def __init__(
        self,
        bucket_path,
        image_path,
        sqlite_path,
        metadata_path,
        plate,
        aggregate_strata_cols,
        gene_name,
        pert_name,
        site,
    ):

        self.sqlite_dir = os.path.join(sqlite_path, plate)
        self.sqlite_file = "sqlite:///{}/{}.sqlite".format(self.sqlite_dir, plate)

        self.barcode_file = os.path.join(metadata_path, "barcode_platemap.csv")
        self.platemap_name = (
            pd.read_csv(self.barcode_file)
            .query("Assay_Plate_Barcode == @plate")
            .Plate_Map_Name.values[0]
        )
        self.platemap_path = os.path.join(
            metadata_path, "platemap", "{}.csv".format(self.platemap_name)
        )

        self.platemap_df = pd.read_csv(self.platemap_path)

        self.channels = ["DNA", "ER", "RNA", "AGP", "Mito"]
        self.channel_colors = {
            "DNA": ["blue", np.array([0, 0, 255], dtype=np.uint8)],
            "ER": ["green", np.array([0, 255, 0], dtype=np.uint8)],
            "RNA": ["yellow", np.array([255, 255, 0], dtype=np.uint8)],
            "AGP": ["orange", np.array([255, 150, 0], dtype=np.uint8)],
            "Mito": ["red", np.array([255, 0, 0], dtype=np.uint8)],
        }

        self.aggregate_strata_cols = aggregate_strata_cols
        self.gene_name = gene_name
        self.pert_name = pert_name
        self.site = site

    def load_image_table(self, get_well=True):
        # Prepare sql file for processing
        self.ap = AggregateProfiles(
            self.sqlite_file, strata=self.aggregate_strata_cols, load_image_data=False
        )

        image_query = "select * from image"
        image_df = pd.read_sql(sql=image_query, con=self.ap.conn)

        image_cols = [
            "Image_Metadata_Plate",
            "Image_Metadata_Well",
            "Image_Metadata_Site",
            "Image_PathName_CellOutlines",
            "Image_PathName_NucleiOutlines",
            "Image_PathName_OrigAGP",
            "Image_PathName_OrigDNA",
            "Image_PathName_OrigER",
            "Image_PathName_OrigMito",
            "Image_PathName_OrigRNA",
            "Image_FileName_CellOutlines",
            "Image_FileName_NucleiOutlines",
            "Image_FileName_OrigAGP",
            "Image_FileName_OrigDNA",
            "Image_FileName_OrigER",
            "Image_FileName_OrigMito",
            "Image_FileName_OrigRNA",
        ]

        self.image_file_core_df = image_df.loc[:, image_cols].merge(
            self.platemap_df,
            left_on="Image_Metadata_Well",
            right_on="well_position",
            how="left",
        )

        if get_well:
            self.find_well()

    def find_well(self, gene_name=None, pert_name=None, site=None):
        if not gene_name:
            gene_name = self.gene_name
        if not pert_name:
            pert_name = self.pert_name
        if not site:
            self.set_site(site)

        self.well = (
            self.platemap_df.query("gene_name == @gene_name")
            .query("pert_name == @pert_name")
            .well_position.values[0]
        )

    def set_site(self, site):
        self.site = site

    def set_well(self, well):
        self.well = well

    def get_image_paths(self):
        self.site_subset_df = self.image_file_core_df.query(
            "Image_Metadata_Well == @self.well"
        ).query("Image_Metadata_Site == @self.site")

        # Create an image channel dictionary
        self.channel_file_dict = {}
        for channel in self.channels:
            file_name_col = "Image_FileName_Orig{}".format(channel)
            path_name_col = "Image_PathName_Orig{}".format(channel)

            file_name = self.site_subset_df.loc[:, file_name_col].values[0]
            path_name = self.site_subset_df.loc[:, path_name_col].values[0]

            path_name = path_name.replace(
                "2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad/2016_04_01_a549_48hr_batch1",
                "2015_07_01_Cell_Health_Vazquez_Cancer_Broad/CRISPR_PILOT_B1",
            )

            file_name = os.path.join(path_name, file_name)
            self.channel_file_dict[channel] = file_name

    def read_images(self):
        # Read Images
        self.image_dict = {}
        for channel_idx in range(0, len(self.channels)):
            channel = self.channels[channel_idx]
            image_key = "ch{}_{}".format(channel_idx + 1, channel)
            self.image_dict[image_key] = skimage.io.imread(
                self.channel_file_dict[channel]
            )

    def crop_and_normalize_images(self, low_prop=0.3, high_prop=0.5):
        # Crop and Normalize Images
        mindim_row = np.rint(self.image_dict["ch1_DNA"].shape[0] * low_prop).astype(
            np.int
        )
        maxdim_row = np.rint(self.image_dict["ch1_DNA"].shape[0] * high_prop).astype(
            np.int
        )

        mindim_col = np.rint(self.image_dict["ch1_DNA"].shape[1] * low_prop).astype(
            np.int
        )
        maxdim_col = np.rint(self.image_dict["ch1_DNA"].shape[1] * high_prop).astype(
            np.int
        )

        self.image_dict_cropped = {}
        for image_key in self.image_dict:
            self.image_dict_cropped[image_key] = normalize_image(
                self.image_dict[image_key][mindim_row:maxdim_row, mindim_col:maxdim_col]
            )

    def colorize_images(self):
        self.image_color_dict = {}
        self.image_color_dict_cropped = {}
        for image_key in self.image_dict:
            img = self.image_dict[image_key]
            img_crop = self.image_dict_cropped[image_key]
            channel = image_key.split("_")[1]
            color, color_array = self.channel_colors[channel]
            self.image_color_dict[image_key] = colorize_image(img, color_array)
            self.image_color_dict_cropped[image_key] = colorize_image(
                img_crop, color_array
            )

    def prep_images(self, low_prop=0.3, high_prop=0.5):
        self.get_image_paths()
        self.read_images()
        self.crop_and_normalize_images(low_prop=low_prop, high_prop=high_prop)
        self.colorize_images()

    def plot_images(self, cropped=True, color=True):
        image_dict = self.get_image_dict(cropped=cropped, color=color)

        fig, ax = plt.subplots(nrows=1, ncols=len(self.channels), figsize=(10, 2))
        for channel_idx in range(0, len(self.channels)):
            channel = self.channels[channel_idx]
            image_key = "ch{}_{}".format(channel_idx + 1, channel)
            channel = image_key.split("_")[1]

            ax[channel_idx].imshow(image_dict[image_key], cmap="gray")
            ax[channel_idx].axis("off")
            ax[channel_idx].set_title(channel)

        plt.tight_layout(pad=0)

    def plot_combined_image(self, cropped=True, color=True, factor=2):
        image_dict = self.get_image_dict(cropped=cropped, color=color)

        combined_image = np.zeros(image_dict["ch1_DNA"].shape)
        for image_key in image_dict:
            combined_image += factor * image_dict[image_key].astype(np.uint16)

        combined_image = normalize(combined_image)

        plt.figure(figsize=(7, 5))
        plt.imshow(combined_image)
        plt.axis("off")
        plt.tight_layout()

    def get_image_dict(self, cropped=True, color=True):
        if cropped:
            if color:
                image_dict = self.image_color_dict_cropped
            else:
                image_dict = self.image_dict_cropped
        else:
            if color:
                image_dict = self.image_color_dict
            else:
                image_dict = self.image_dict
        return image_dict


def normalize_image(img):
    """
    Function copied from https://github.com/jr0th/segmentation/blob/master/visualization/CellArt.ipynb
    """
    # normalize to [0,1]
    percentile = 99
    high = np.percentile(img, percentile)
    low = np.percentile(img, 100 - percentile)

    img = np.minimum(high, img)
    img = np.maximum(low, img)

    # gives float64, thus cast to 8 bit later
    img = (img - low) / (high - low)

    img = skimage.img_as_ubyte(img)
    return img


def colorize_image(img, col):
    """
    Function copied from https://github.com/jr0th/segmentation/blob/master/visualization/CellArt.ipynb
    """
    # rescale image
    img_float = img.astype(np.float)
    img_float = img_float / 255

    # colorize
    img_col_float = np.reshape(img_float, img_float.shape + (1,)) * col
    img_col_byte = img_col_float.astype(np.uint8)

    return img_col_byte
