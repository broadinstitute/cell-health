#!/bin/bash
#
# Gregory Way 2019
# Predicting Cell Health
#
# Generate Audits - Assess Replication Correlation

#############################
# Step 1 - Setup Constants
#############################
PLATES=../../scratch/processed_plates.txt
PLATE_MAPS=../../scratch/plate_map_names.txt
BATCH_ID=CRISPR_PILOT_B1
SAMPLE_PLATE_ID=SQ00014610
PROJECT_NAME=2015_07_01_Cell_Health_Vazquez_Cancer_Broad

# Create directory
output_dir=../../scratch/${BATCH_ID}/audit/
mkdir -p output_dir

##############################
# Step 2 - Audit Replicate Correlations
##############################
# Here, we care about the replicate correlation of:
# 1. Controls
# 2. Treatment
# 3. CRISPR guides targetting the same genes

data_abbrev=( 'control' 'treatment' 'crispr')
for data_type in "${data_abbrev[@]}"
do
  if [ $data_type == "treatment" ]; then
    subset_type="trt"
    output_string=""
    group_by_variables=Metadata_Plate_Map_Name,Metadata_Well,Metadata_gene_name,Metadata_pert_name,Metadata_broad_sample,Metadata_cell_line,Metadata_pert_id,Metadata_pert_mfc_id,Metadata_cell_id,Metadata_broad_sample_type,Metadata_pert_type
  elif [ $data_type == "control" ]; then
    subset_type="control"
    output_string="_control"
    group_by_variables=Metadata_cell_line,Metadata_Well,Metadata_pert_name
  else
    subset_type="trt"
    output_string="_guide"
    group_by_variables=Metadata_cell_line,Metadata_Plate_Map_Name,Metadata_gene_name,Metadata_pert_type
  fi
  parallel \
    --no-run-if-empty \
    --eta \
    --joblog ../../log/${BATCH_ID}/audit_control.log \
    --results ../../log/${BATCH_ID}/audit_control \
    --keep-order \
    ./audit.R \
    --batch_id ${BATCH_ID} \
    --plate_map_name {1} \
    --suffix _normalized_variable_selected.csv \
    --subset \"Metadata_broad_sample_type == \'\'\'${subset_type}\'\'\'\" \
    --output ${output_dir}/{1}_audit${output_string}.csv \
    --output_detailed ${output_dir}/{1}_audit_detailed${output_string}.csv \
    --group_by ${group_by_variables} :::: ${PLATE_MAPS}
done
