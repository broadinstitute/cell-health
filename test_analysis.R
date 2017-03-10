library(tidyverse)
df <- read_csv("../../scratch/CRISPR_PILOT_B1/CRISPR_PILOT_B1_collapsed.csv")


df %>% filter(Metadata_pert_name != "EMPTY") %>% write_csv("../../scratch/CRISPR_PILOT_B1/CRISPR_PILOT_B1_collapsed_no_empty.csv")