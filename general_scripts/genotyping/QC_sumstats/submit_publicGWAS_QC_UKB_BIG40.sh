#!/bin/bash
# Prepare public GWAS sumstats from UKB to be used as base datasets for PRS analyses or for genetic correlation analyses (LDSC)
#---------------------------------------------------------------------- 
# From UKB BIG40 study
## https://open.win.ox.ac.uk/ukbiobank/big40/
#---------------------------------------------------------------------- 
primary_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/downloaded_data/UKB/BIG40/
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/QC/UKB_BIG40/
analysis_dir=/export/home/acarrion/acarrion/projects/general_scripts/genotyping/QC_sumstats/
#----------------------------------------------------------------------
mkdir -p ${working_dir}

# Select regions to format: 
# Rscript --verbose ${primary_dir}select_phenos2download.R

cd ${primary_dir}
# -e corpus -e operc -e planum_temp 
files2format=$(awk '{print $1}' IDPs_downloaded.txt | sed 's/"//g' | grep -v Pheno | sort | uniq)

#
cd ${analysis_dir}
mkdir -p sge_jobs
cd sge_jobs
for base_pheno in ${files2format}
 do
  base_file=${base_pheno}.txt.gz
  # extract phenotype name from file
  if [ ! -f ${working_dir}${base_pheno}.QC.txt.gz ]
  then
    echo Submit ${base_pheno}
    qsub ${analysis_dir}prepare_publicGWAS_QC_UKB_BIG40.sh ${base_file}
  else
    echo Skip, ${base_pheno} already QCd.
  fi
 done
