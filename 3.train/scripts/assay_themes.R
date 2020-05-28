# Specific themes, colors, and text to be held consistent in visualizing assay model performance

library(ggplot2)

measurement_levels <- c(
    "cell_viability",
    "dna_damage",
    "g1_phase",
    "s_phase",
    "g2_phase",
    "early_mitosis",
    "mitosis",
    "late_mitosis",
    "cell_cycle_count",
    "shape",
    "apoptosis",
    "death",
    "ros",
    "other"
)

measurement_colors <- c(
    "cell_viability"=  "#b2df8a",
    "dna_damage" = "#fb9a99",
    "g1_phase" = "#fdbf6f",
    "s_phase" = "#cab2d6",
    "g2_phase" = "#ff7f00",
    "early_mitosis" = "#ccece6",
    "mitosis" = "#41ae76",
    "late_mitosis" = "#005824",
    "cell_cycle_count" = "#1f78b4",
    "shape" = "#6a3d9a",
    "apoptosis" = "#a6cee3",
    "death" = "#33a02c",
    "ros" = "#e31a1c",
    "other" = "black"
)

measurement_labels <- c(
    "cell_viability"=  "Cell Viability",
    "dna_damage" = "DNA Damage",
    "g1_phase" = "G1 Phase",
    "s_phase" = "S Phase",
    "g2_phase" = "G2 Phase",
    "early_mitosis" = "Early Mitosis",
    "mitosis" = "Mitosis",
    "late_mitosis" = "Late Mitosis",
    "cell_cycle_count" = "Cell Cycle Count",
    "shape" = "Shape",
    "apoptosis" = "Apoptosis",
    "death" = "Death",
    "ros" = "Reactive Oxygen Species",
    "other" = "CRISPR Efficiency"
)

dye_levels <- c(
    "hoechst",
    "edu",
    "hoechst_edu",
    "hoechst_edu_ph3",
    "hoechst_gh2ax",
    "hoechst_edu_gh2ax",
    "hoechst_edu_ph3_gh2ax",
    "draq",
    "draq_caspase",
    'cell_rox',
    "dpc",
    "qc"
)

dye_colors <- c(
    "hoechst" = "#639B94",
    "edu" = "#E45242",
    "hoechst_edu" = "#73414b",
    "hoechst_edu_ph3" = "#7B9C32",
    "hoechst_gh2ax" = "#535f52",
    "hoechst_edu_gh2ax" = "#e37a48",
    "hoechst_edu_ph3_gh2ax" = "#E2C552",
    "draq" = "#FF6699",
    "draq_caspase" = "#7f4a72",
    'cell_rox' = "#E9DFC3",
    "dpc" = "#CA662D",
    "qc" = "black"
)

dye_labels <- c(
  "hoechst" = "Hoechst",
  "edu" = "EdU",
  "hoechst_edu" = "Hoechst + EdU",
  "hoechst_edu_ph3" = "Hoechst + EdU + PH3",
  "hoechst_gh2ax" = "Hoechst + gH2AX",
  "hoechst_edu_gh2ax" = "Hoechst + EdU + gH2AX",
  "hoechst_edu_ph3_gh2ax" = "Hoechst + EdU + PH3 + gH2AX",
  "draq" = "DRAQ7",
  "draq_caspase" = "DRAQ7 + Caspase",
  'cell_rox' = "CellROX",
  "dpc" = "DPC (Shape)",
  "qc" = "CRISPR Efficiency"
)

cell_line_labels <- c(
  "A549" = "A549",
  "ES2" = "ES2",
  "HCC44" = "HCC44"
)

cell_line_colors <- c(
  "A549" = "#861613",
  "ES2" = "#1CADA8",
  "HCC44" = "#2A364D"
)

dye_theme <- theme(
    axis.title = element_text(size = 4),
    axis.title.x = element_text(margin = margin(0, 0, 0, 0)),
    axis.title.y = element_text(margin = margin(0, 0, 0, 0)),
    axis.text = element_text(size = 4),
    axis.ticks = element_line(size = 0.1),
    axis.ticks.length = unit(0.05, "cm"),
    plot.title = element_text(size = 5, margin = margin(0, 0, 0, 0)),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"),
    legend.title = element_text(size = 3.5),
    legend.text = element_text(size = 3.5),
    legend.key.height = unit(0.2, "line"),
    legend.key.width = unit(-0.2, "line"),
    panel.grid.minor = element_line(size = 0.1),
    panel.grid.major = element_line(size = 0.2),
    legend.margin = margin(-5, 7, 0, 0),
    legend.box.margin=margin(0, -5, 5, -5)
)
