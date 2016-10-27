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

for (cell_line in c("A549", "ES2", "HCC44")) {
  platemap %>%
    mutate(cell_line = cell_line) %>%
    write_csv(sprintf("../../metadata/CRISPR_PILOT_B1/platemap/DEPENDENCIES1_%s.csv", cell_line))
  
  platemap %>%
    mutate(cell_line = cell_line) %>%
    write_tsv(sprintf("../../metadata/CRISPR_PILOT_B1/platemap/DEPENDENCIES1_%s.txt", cell_line))
}

