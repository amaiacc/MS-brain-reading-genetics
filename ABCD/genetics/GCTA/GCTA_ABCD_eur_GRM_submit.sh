#!/bin/bash

# Build GRM to compute heritability estimates using GCTA: ABCD sample
# date: 27.11.2020
# by: Amaia Carrion Castillo

# using gcta version 1.92beta
# http://cnsgenomics.com/software/gcta/#MakingaGRM

# first: clean input data

#----------------------------------------------------------------------
# Define paths and directories
#----------------------------------------------------------------------
resource_dir=/export/home/acarrion/acarrion/projects/resources/bin/
gcta_dir=${resource_dir}gcta_1.93.2beta/
# working directory
project_dir=/export/home/acarrion/acarrion/projects/resources/datasets/ABCD/
input_dir=${project_dir}data/working_data/genotyping/preImputation_QC/
working_dir=${project_dir}data/working_data/genotyping/GCTA/grm/
mkdir -p ${working_dir}


scripts_dir=${project_dir}/scripts/genetics/GCTA/
cd ${scripts_dir}/sge_jobs/
# build grm
for root in ABCD_release_3.0_QCed_EUR.EUR6sd_unrelated.cleaned ABCD_release_3.0_QCed_EUR.EUR3sd_unrelated.cleaned
do
 qsub ${scripts_dir}GCTA_GRM_allchrs.sh ${root}
done
