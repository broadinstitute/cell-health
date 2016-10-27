library(tidyverse)
library(stringr)

platemap <- 
  read_csv("../../metadata/CRISPR_PILOT_B1/platemap/DEPENDENCIES1_raw.csv") %>%
  gather(WellCol, pert_name, -Row) %>%
  rename(WellRow = Row) %>%
  rowwise() %>%
  mutate(well_position = sprintf("%s%02d", WellRow, as.integer(WellCol))) %>%
  mutate(gene_name = str_split(pert_name, "-")[[1]][1]) %>%
  ungroup() %>%
  select(WellRow, WellCol, well_position, gene_name, pert_name) %>%
  arrange(well_position) 

platemap %>%
  write_csv("../../metadata/CRISPR_PILOT_B1/platemap/DEPENDENCIES1.csv")

platemap %>%
  write_tsv("../../metadata/CRISPR_PILOT_B1/platemap/DEPENDENCIES1.txt")