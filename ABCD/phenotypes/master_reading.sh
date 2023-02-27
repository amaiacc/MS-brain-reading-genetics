#!/bin/bash

module load R/4.0.3
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/ABCD/data/working_data/phenotypes/
scripts_dir=/export/home/acarrion/acarrion/projects/resources/datasets/ABCD/scripts/phenotypes/

cd ${scripts_dir}
#dos2unix *

config_file=${scripts_dir}reading.config

# get all keys from config file --> to bash
keys=( $(grep -oP '\w+(?==)' reading.config) )
for (( i=0; i < ${#keys[@]}; i++ )); do
    printf "%d\t%s\n" $i "${keys[i]}"
done
echo
source ${config_file}
for var in "${keys[@]}"; do
    printf "%s\t=> %s\n" "$var" "${!var}"
done

mkdir -p ${scripts_dir}/sge
cd ${scripts_dir}/sge

#--------------------------------
# 0 baseline models
## Rscript 0_baseline.R ${config_file}
#--------------------------------

#--------------------------------
# 1 ROIs
#--------------------------------
# baseline model for reading, i.e. make sure model diagnostics look fine
## Rscript 0_rois_baseline.R ${config_file}
#--------------------------------
# reading ~ ROIs

## parallelize per model
# - arguments are: config file, test model (row)

models=${working_dir}/${v}/${s}/${pheno}/baseline_${atlas}/trim${trim_val}/tables/baseline_models_table.csv
n=$(wc $models -l | awk '{print $1}') # run in parallel
n2=$(($n-1))


for ((i=1;i<=n2;i++)); do
 echo ${i}
 # run if output not present
 qsub ${scripts_dir}models_submit.sh ${config_file} ${i}

done


# clean sge job files: remove if no analysis run
for f in $(ls *models.e*); do n=$(grep REML $f | wc -l ); echo $f $n; done > tmp
files2rm=$(grep " 0$" tmp | awk '{print $1}' | sed 's/abcd_models.e//g')
for f in $files2rm; do rm *$f; done
rm tmp
#--------------------------------

#--------------------------------
# 2 ROIs - parse output
#--------------------------------
cd ${scripts_dir}
Rscript 2_rois_models_parse_main.R

#--------------------------------
# 3 ROIs - visualization
#--------------------------------



#--------------------------------
# 1 follow up on RH rois2follow
#--------------------------------

## parallelize per model
# - arguments are: config file, test model (row)

models=${working_dir}/${v}/${s}/${pheno}/baseline_${atlas}/trim${trim_val}/tables/baseline_models_table.csv
n=$(cat $models | grep -v -e picvocab -e pea_wiscv -e fluidcomp | wc -l | awk '{print $1}') # run in parallel
n2=$(($n-1))

for ((i=1;i<=n2;i++)); do
 echo ${i}
 # run if output not present
 qsub ${scripts_dir}models_submit_RH.sh ${config_file} ${i}

done
