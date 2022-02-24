#!/bin/bash

#----------------------------------------------------------------------
# Run pre-imputation QC for ABCD dataset
#----------------------------------------------------------------------
# previous: 
# - plink files downloaded from 
# NOTE: zCall recalls all NC's --> zCall output should not be used for QC of samples and SNPs
# 1- define bad quality samples and SNPs using plink (from GenomeStudio output without recalling -> converted to plink format in GenomeStudio_zCall_rareSNPrecalling.sh)
# 2- remove bad quality samples and SNPs in plink (--remove and --exclude)
#----------------------------------------------------------------------
module load plink/6.15
#
qc_scripts=/export/home/acarrion/acarrion/projects/general_scripts/genotyping/QC/
#
primary_dir=/export/home/acarrion/acarrion/projects/resources/datasets/ABCD/data/primary_data/ABCDgenomicsDataV30/genomics_sample03/ABCD_genotype/genotype_QCed/

# working directory
project_dir=/export/home/acarrion/acarrion/projects/resources/datasets/ABCD/
geno_scripts=${project_dir}scripts/genetics/QC/
working_dir=${project_dir}data/working_data/genotyping/preImputation_QC/
# create working dir
mkdir -p ${working_dir} ${geno_scripts}
# copy latest scripts to geno_scripts
cp ${qc_scripts}/*R ${geno_scripts}
cp ${qc_scripts}/*sh ${geno_scripts}


#----------------------------------------------------------------------
# run QC per genotyping batch
cd ${project_dir}scripts
root=ABCD
#batch=release_2.0.1_r1
batch=release_3.0_QCed
bash ${geno_scripts}/preImputation_batchQC_ABCD.sh ${project_dir} ${root} ${batch} ${primary_dir}
#----------------------------------------------------------------------
