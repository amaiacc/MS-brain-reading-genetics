#!/bin/bash
# PRS
#$ -N PRS_parse_dataset
#$ -cwd
#$ -q short.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m beas
#----------------------------------------------------------------------
# Run PRS analysis on the ABCD dataset, using public GWAS sumstats for cognitive measures as base datasets
#----------------------------------------------------------------------
bin_dir=/export/home/acarrion/acarrion/projects/resources/bin/
prsice_dir=${bin_dir}PRSice_v2.2.12/
PATH=${prisce_dir}:${PATH}
#----------------------------------------------------------------------
# Load modules, define directories, root files, and PATH
#----------------------------------------------------------------------
module load plink/6.15
module load R/4.0.3
#----------------------------------------------------------------------
# define variables
root=ABCD
batch=release_3.0_QCed
data=${root}_${batch}
base_project=cognitive
#----------------------------------------------------------------------
# build parameters
#----------------------------------------------------------------------
add_dir=/resources/datasets/
analysis_dir=/export/home/acarrion/acarrion/projects/${add_dir}${root}/scripts/genetics/PRS/

base_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/QC/${base_project}/

project_dir=/export/home/acarrion/acarrion/projects/${add_dir}${root}/data/working_data/genotyping/
target_dir=${project_dir}/imputation/${data}/MichiganImputationServer/imputed_genotypes/clean/
if [ ! -d ${target_dir} ]
 then
 target_dir=${project_dir}imputed_genotypes/clean/
fi
working_dir=${project_dir}/prsice/

mkdir -p ${working_dir}/${base_project}/
#----------------------------------------------------------------------
cd ${analysis_dir}/sge


#Rscript ${analysis_dir}PGS_QC.R ${analysis_dir}base_${base_project}.config --verbose
Rscript ${analysis_dir}PGS_lmm.R ${analysis_dir}base_${base_project}.config --verbose