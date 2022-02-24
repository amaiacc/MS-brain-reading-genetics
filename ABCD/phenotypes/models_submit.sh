#!/bin/bash
# Run lmer on ABCD reading ~ sMRI ROIs
#$ -N abcd_models
#$ -cwd
#$ -q short.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m as
#----------------------------------------------------------------------
module load R/4.0.3

#subset=${1}
#trim_val=${2}
#pheno=${3}
#atlas=${4}
config_file=${1}
i=${2}


cd /export/home/acarrion/acarrion/projects/resources/datasets/ABCD/scripts/phenotypes/

#Rscript --verbose 1reading_rois_models.R  ${subset} ${trim_val} ${pheno} ${atlas} ${i}
Rscript --verbose 1_rois_models.R ${config_file} ${i}