#!/bin/bash
# Prepare public GWAS sumstats from UKB to be used as base datasets for PRS analyses or for genetic correlation analyses (LDSC)
#---------------------------------------------------------------------- 
# From UKB BIG40 study
## https://open.win.ox.ac.uk/ukbiobank/big40/
#---------------------------------------------------------------------- 
primary_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/ldsc/ieu-openGWAS/ # mungedSumstats, original files are in vcf format
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/QC/ieu-openGWAS/
#----------------------------------------------------------------------
# define variables
root=ABCD
batch=release_3.0_QCed
data=${root}_${batch}
base_project=sibGWAS_Howe2022
#----------------------------------------------------------------------
# build parameters
#----------------------------------------------------------------------
add_dir=/resources/datasets/
analysis_dir=/export/home/acarrion/acarrion/projects/${add_dir}${root}/scripts/genetics/${base_project}/
#----------------------------------------------------------------------


#----------------------------------------------------------------------
mkdir -p ${working_dir}

# Select regions to format: 
# Rscript --verbose ${primary_dir}select_phenos2download.R

cd ${primary_dir}
# -e corpus -e operc -e planum_temp 
files2format=$(ls *tsv.gz)

#
cd ${analysis_dir}
for base_file in ${files2format}
 do
  base_pheno=$(echo ${base_file} | sed 's/.tsv.gz//g')
  # extract phenotype name from file
  if [ ! -f ${working_dir}${base_pheno}.QC.txt.gz ]
  then
    echo Submit ${base_pheno}
    qsub ${analysis_dir}prepare_${base_project}_QC.sh ${base_file}
  else
    echo Skip, ${base_pheno} already QCd.
  fi
 done
