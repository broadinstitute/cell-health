# Predicting Cell Health
# Tim Becker and Gregory Way
#
# The following scripts were performed to generate profiles for the Cell Health Project

#############################
# Step 1 - Setup Constants
#############################
PLATES=../../scratch/processed_plates.txt
PLATE_MAPS=../../scratch/plate_map_names.txt
BATCH_ID=CRISPR_PILOT_B1
SAMPLE_PLATE_ID=SQ00014610
PROJECT_NAME=2015_07_01_Cell_Health_Vazquez_Cancer_Broad
BROAD_NFS=/imaging/analysis
AWS_NFS=/home/ubuntu/bucket/projects/
BUCKET=imaging-platform
MAXPROCS=3

#############################
# Step 2 - Setup Directory Structure
#############################
mkdir -p ${AWS_NFS}/${PROJECT_NAME}/workspace/
cd ${AWS_NFS}/${PROJECT_NAME}/workspace/
mkdir -p log/${BATCH_ID}

# For use with cytominer-database
mkdir ${AWS_NFS}/ebs_tmp

#############################
# Step 3 - Add GitHub Software
#############################
cd ${AWS_NFS}/${PROJECT_NAME}/workspace/
mkdir software && cd software

git clone git@github.com:broadinstitute/cytominer_scripts.git
git clone git@github.com:broadinstitute/$PROJECT_NAME

##############################
# Step 4 - Load Data
##############################
cd ${AWS_NFS}/${PROJECT_NAME}/workspace

# Download all data (except for sqlite database) from backend
mkdir backend && cd backend
aws s3 sync --exclude "*.sqlite" "s3://imaging-platform/projects/2015_07_01_Cell_Health_Vazquez_Cancer_Broad/workspace/backend/" .
cd ..

# Download the remaining folders and data
DIRS=( collated metadata results log scratch audit parameters )
for DIR in "${DIRS[@]}"
do
  mkdir $DIR
  cd $DIR
  aws s3 sync "s3://imaging-platform/projects/2015_07_01_Cell_Health_Vazquez_Cancer_Broad/workspace/$DIR" . && cd ..
done

##############################
# Step 5 - Generate Plate Maps
##############################
cd ${AWS_NFS}/${PROJECT_NAME}/workspace/software/cytominer_scripts/

# generate PLATE_MAPS
csvcut -c Plate_Map_Name \
  ../../metadata/${BATCH_ID}/barcode_platemap.csv | \
  tail -n +2|sort|uniq > \
  ${PLATE_MAPS}

# generate PLATES
find ../../backend/ -type f -name "*.csv"|cut -d"/" -f5|sort|uniq > ${PLATES}

##############################
# Step 5 - Process Plates
##############################
# process all plates
mkdir -p  ../../log/${BATCH_ID}/
parallel \
  --max-procs ${MAXPROCS} \
  --eta \
  --joblog ../../log/${BATCH_ID}/collate.log \
  --results ../../log/${BATCH_ID}/collate \
  --files \
  --keep-order \
  ./collate.R \
  -b ${BATCH_ID} \
  --plate {1} \
  -c ingest_config.ini \
  --tmpdir ~/ebs_tmp \
  -d \
  -r s3://${BUCKET}/projects/${PROJECT_NAME}/workspace :::: ${PLATES}

parallel ./process.sh \
  -b ${BATCH_ID} \
  -p {1} \
  -t ~/ebs_tmp/ :::: ${PLATES}

# annotate
parallel ./annotate.R \
  -b ${BATCH_ID} \
  -p {1} \
  -d \
  -m genetic :::: ${PLATES}

# normalize
# NOTE: don't escape quotes if not using parallel
parallel ./normalize.R \
  -b ${BATCH_ID} \
  -p {1} \
  -s \"Metadata_broad_sample_type == \'\'\'control\'\'\'\" :::: ${PLATES}

# create samples to do feature selection
mkdir -p ../../parameters/${BATCH_ID}/sample/

# sample normalized data
./sample.R \
  -b ${BATCH_ID} \
  -f "_normalized.csv$" \
  -n 2 \
  -o ../../parameters/${BATCH_ID}/sample/${BATCH_ID}_normalized_sample.feather

# sample unnormalized data
./sample.R \
  -b ${BATCH_ID} \
  -f "_augmented.csv$" \
  -n 2 \
  -o ../../parameters/${BATCH_ID}/sample/${BATCH_ID}_augmented_sample.feather

# replicate_correlation
./preselect.R \
  -b ${BATCH_ID} \
  -i ../../parameters/${BATCH_ID}/sample/${BATCH_ID}_normalized_sample.feather \
  -r replicate_correlation \
  -s "Metadata_broad_sample_type == '''trt'''" \
  -n 2

# correlation_threshold
./preselect.R \
  -b ${BATCH_ID} \
  -i ../../parameters/${BATCH_ID}/sample/${BATCH_ID}_normalized_sample.feather \
  -r correlation_threshold

# variance_threshold
./preselect.R \
  -b ${BATCH_ID} \
  -i ../../parameters/${BATCH_ID}/sample/${BATCH_ID}_augmented_sample.feather \
  -r variance_threshold \
  -s "Metadata_broad_sample_type == '''control'''"

# manually remove some features
echo "variable" > ../../parameters/${BATCH_ID}/variable_selection/manual.txt

head -1 \
  ../../backend/${BATCH_ID}/${SAMPLE_PLATE_ID}/${SAMPLE_PLATE_ID}.csv \
  |tr "," "\n"|grep -v Meta|grep -E -v 'Granularity_14|Granularity_15|Granularity_16|Manders|RWC' >> \
  ../../parameters/${BATCH_ID}/variable_selection/manual.txt

# feature select
parallel ./select.R \
  -b ${BATCH_ID} \
  -p {1} \
  -r variance_threshold,replicate_correlation,correlation_threshold,manual :::: ${PLATES}

# create gct for each plate
parallel ./csv2gct.R \
  ../../backend/${BATCH_ID}/{1}/{1}_normalized_variable_selected.csv \
  -o ../../backend/${BATCH_ID}/{1}/{1}_normalized_variable_selected.gct :::: ${PLATES}

# audit each plate_map for replicate reproducibility
mkdir -p ../../scratch/${BATCH_ID}/audit/

# only treatments
parallel ./audit.R \
  -b ${BATCH_ID} \
  -m {1} \
  -s \"Metadata_broad_sample_type == \'\'\'trt\'\'\'\" \
  -o ../../scratch/${BATCH_ID}/audit/{1}_audit.csv \
  -p Metadata_pert_name :::: ${PLATE_MAPS}

# only controls
parallel ./audit.R \
  -b ${BATCH_ID} \
  -m {1} \
  -s \"Metadata_broad_sample_type == \'\'\'control\'\'\'\" \
  -o ../../scratch/${BATCH_ID}/audit/{1}_audit_control.csv \
  -p Metadata_pert_name :::: ${PLATE_MAPS}

# only treatments
parallel \
  --no-run-if-empty \
  --eta \
  --joblog ../../log/${BATCH_ID}/audit.log \
  --results ../../log/${BATCH_ID}/audit \
  --files \
  --keep-order \
  ./audit.R \
  -b ${BATCH_ID} \
  -m {1} \
  -f _normalized_variable_selected.csv \
  -s \"Metadata_broad_sample_type == \'\'\'trt\'\'\'\" \
  -o ../../scratch/${BATCH_ID}/audit/{1}_audit.csv \
  -l ../../scratch/${BATCH_ID}/audit/{1}_audit_detailed.csv \
  -p Metadata_Plate_Map_Name,Metadata_Well,Metadata_gene_name,Metadata_pert_name,Metadata_broad_sample,Metadata_cell_line,Metadata_pert_id,Metadata_pert_mfc_id,Metadata_cell_id,Metadata_broad_sample_type,Metadata_pert_type :::: ${PLATE_MAPS}
`
# only controls
parallel \
  --no-run-if-empty \
  --eta \
  --joblog ../../log/${BATCH_ID}/audit_control.log \
  --results ../../log/${BATCH_ID}/audit_control \
  --keep-order \
  ./audit.R \
  -b ${BATCH_ID} \
  -m {1} \
  -f _normalized_variable_selected.csv \
  -s \"Metadata_broad_sample_type == \'\'\'control\'\'\'\" \
  -o ../../scratch/${BATCH_ID}/audit/{1}_audit_control.csv \
  -p Metadata_Well :::: ${PLATE_MAPS}

# collapse
mkdir -p ../../scratch/${BATCH_ID}/collapsed

parallel ./collapse.R \
  -b ${BATCH_ID} \
  -m {1} \
  -o ../../scratch/${BATCH_ID}/collapsed/{1}_collapsed.csv :::: ${PLATE_MAPS}

csvstack \
  ../../scratch/${BATCH_ID}/collapsed/*_collapsed.csv > \
  ../../scratch/${BATCH_ID}/${BATCH_ID}_collapsed.csv

# convert collapsed to gct
./csv2gct.R \
  ../../scratch/${BATCH_ID}/${BATCH_ID}_collapsed.csv \
  -o ../../scratch/${BATCH_ID}/${BATCH_ID}_collapsed.gct

########################################
# Generate Heterogeneity Features
########################################
parallel \
  --no-run-if-empty \
  --eta \
  --joblog ../../log/${BATCH_ID}/audit_control.log \
  --results ../../log/${BATCH_ID}/audit_control \
  --keep-order \
  ./capture_heterogeneity.R \
  --sqlite_file ../../backend/${BATCH_ID}/{1}/{1}.sqlite \
  --num_comp 3000 {1} \
  --output ../../backend/${BATCH_ID}/{1}/{1}_heterogeneity_feature.csv :::: ${PLATE_MAPS}

# sync up to S3 from EC2 instance
# on ec2
aws s3 sync \
  --exclude "*" \
  --include "*.gct" \
  ../../backend/${BATCH_ID}/ \
  s3://imaging-platform/projects/${PROJECT_NAME}/workspace/backend/${BATCH_ID}/

# sync down from S3 on Broad instance
# on prefix
aws s3 sync \
  --exclude "*" \
  --include "*.gct" \
  s3://imaging-platform/projects/${PROJECT_NAME}/workspace/backend/${BATCH_ID}/ \
  ${BROAD_NFS}/${PROJECT_NAME}/workspace/backend/${BATCH_ID}/

# compare plates
# on prefix
cd ${BROAD_NFS}/${PROJECT_NAME}/workspace/software/cytominer_scripts/

mkdir -p ../../scratch/${BATCH_ID}/compare_plates

cat ${PLATE_MAPS} | \
  xargs -P 28 -I % \
  ./compare_plates.R \
  -b ${BATCH_ID} \
  -m % \
  -o ../../scratch/${BATCH_ID}/compare_plates

# sync up to S3 from EC2 instance
# on ec2
aws s3 sync \
  --exclude "*" \
  --include "*.gct" \
  ../../scratch/ \
  s3://imaging-platform/projects/${PROJECT_NAME}/workspace/scratch/

# sync down from S3 on Broad instance
# on prefix
aws s3 sync \
  --exclude "*" \
  --include "*.gct" \
  s3://imaging-platform/projects/${PROJECT_NAME}/workspace/scratch/ \
  ${BROAD_NFS}/${PROJECT_NAME}/workspace/scratch/

# sig_introspect
# on prefix
cd ${BROAD_NFS}/${PROJECT_NAME}/workspace/software/cytominer_scripts/

mkdir -p ../../scratch/${BATCH_ID}/sig_introspect

rum -q local sig_introspect_tool \
  --sig_score ../../scratch/${BATCH_ID}/${BATCH_ID}_collapsed.gct \
  -o ../../scratch/${BATCH_ID}/sig_introspect

# pcl analysis

# create moa sets

cd ${BROAD_NFS}/${PROJECT_NAME}/workspace/software/${PROJECT_NAME}

./create_perturbagen_set.R

# UPDATE the path
INTROSPECT_OUTPUT_DIR_TAG=jan04/my_analysis.mortar.tools.SigIntrospectTool.2017010413260591/
INTROSPECT_OUTPUT_DIR=${BROAD_NFS}/${PROJECT_NAME}/workspace/scratch/${BATCH_ID}/sig_introspect/${INTROSPECT_OUTPUT_DIR_TAG}

cd ${INTROSPECT_OUTPUT_DIR}
mkdir for_pcl
cp ps_n10752x10752.gctx for_pcl/ps_A549_n10752x10752.gctx
cd for_pcl/

# aggregate across dose - pick more extreme of 33rd and 67th percentile
# UPDATE the file names
echo \
"ds = parse_gctx('ps_A549_n10752x10752.gctx');
agg_fun = @(x, dim) max_quantile(x, 33, 67, dim);
ds_agg = ds_aggregate(ds,...
         'row_fields', {'pert_id'},...
         'col_fields', {'pert_id'},...
         'fun',  agg_fun);

mkgctx('ps_A549.gctx', ds_agg)" > agg.m

matlab -nosplash -nodesktop -nojvm < agg.m

# create row_meta and col_meta in order to run pcl_eval
# UPDATE the file names
echo \
"library(magrittr)
library(dplyr)
df <- roller::parse.gctx('ps_A549_n1571x1571.gctx')

format_id <- function(desc) {
    stopifnot('id' %in% names(desc))
    desc %>% select_(.dots = c('id', setdiff(names(desc), c('id'))))
}
roller::write.tbl(df@rdesc %>% format_id() %>% rename(rid = id), sprintf('row_meta_n%d.txt', nrow(df@rdesc)))
roller::write.tbl(df@cdesc %>% format_id() %>% rename(cid = id), sprintf('col_meta_n%d.txt', nrow(df@cdesc)))
" > write_meta.R

Rscript -e "source('write_meta.R')"

mkdir sandbox

# UPDATE the file names
mv ps_A549_n1571x1571.gctx row_meta_n1571.txt col_meta_n1571.txt sandbox/

mkdir -p ${BROAD_NFS}/${PROJECT_NAME}/workspace/scratch/${BATCH_ID}/pcleval_all

# run pcl_eval
# UPDATE the file names
rum -q local sig_pcleval_tool \
   --inpath ${INTROSPECT_OUTPUT_DIR}/for_pcl/sandbox \
   --pcl    ${BROAD_NFS}/${PROJECT_NAME}/workspace/scratch/${BATCH_ID}/perturbagen_sets/moa_set_all.gmt \
   -o       ${BROAD_NFS}/${PROJECT_NAME}/workspace/scratch/${BATCH_ID}/pcleval_all
