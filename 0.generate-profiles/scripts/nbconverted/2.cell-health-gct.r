
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))

# Load function
gct_file <- file.path("..", "..", "cytominer_scripts", "write_gct.R")
source(gct_file)

# Load Cell Health Labels and Recode Metadata
label_file <- file.path("data", "labels", "normalized_cell_health_labels.tsv")
labels_df <- readr::read_tsv(label_file, col_types=readr::cols()) %>%
    dplyr::rename("Metadata_cell_id" = "cell_id",
                  "Metadata_guide" = "guide",
                  "Metadata_plate_name" = "plate_name",
                  "Metadata_well_col" = "well_col",
                  "Metadata_well_row" = "well_row")

head(labels_df, 3)

channels <- NULL
create_row_annotations <- TRUE
feature_regex <- "^cc_|vb_"
output <- file.path("data", "labels", "normalized_cell_health_labels.gct")

write_gct(x = labels_df,
          path = output,
          channels = channels,
          create_row_annotations = create_row_annotations,
          feature_regex = feature_regex)
