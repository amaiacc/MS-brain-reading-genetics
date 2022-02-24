#!/bin/bash
# PRS
#$ -N ldsc_format_UKB
#$ -cwd
#$ -q short.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m as
#--------
# format summary statistics create input files for LDSC
#--------

#--------
# get arguments
#--------
p=${1}
pheno_file=${2}
#--------

#--------
# load modules and define tools/directories with resources
module load python/python2.7
PATH=/export/home/acarrion/acarrion/projects/resources/github/ldsc/:${PATH}
resource_dir=/export/home/acarrion/acarrion/projects/resources/
## info on LD score formats/ estimation
# https://github.com/bulik/ldsc/wiki/LD-File-Formats
# https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial
ldscores_dir=${resource_dir}/LDscores/
#--------
# project specific parameters, point to data location and working directory
assoc_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/QC/UKB_BIG40/
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/ldsc/UKB_BIG40/
mkdir -p ${working_dir}
cd ${working_dir}
#--------

#--------
## sumstats format: https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format
# required columns:
##    SNP -- SNP identifier (e.g., rs number)
##  N -- sample size (which may vary from SNP to SNP).
##  Z -- z-score. Sign with respect to A1 (warning, possible gotcha)
##  A1 -- first allele (effect allele) 
##  A2 -- second allele (other allele)
# Note that ldsc filters out all variants that are not SNPs and strand-ambiguous SNPs.
#--------
# GWAS-es from the UKB were run using BGENIE
# BGENIE: 'the regression model we code the first and second alleles as 0 and 1 respectively, so the beta coefficient refers to the effect of having an extra copy of the second allele.'
# a_2 The character code for the first allele (a string). --> effect allele (because BIG40 released betas as direction of a_2!)
# a_1 The character code for the second allele (a string). --> not-effect allele
## https://jmarchini.org/bgenie-usage/
#--------
##  QCd sumstats current format:
#chr rsid pos a_1 a_0 beta se pval(-log10) P
#--------



p_file=${assoc_dir}${p}.QC.txt.gz

if [ ! -f ${p}.sumstats.gz ]
 then
  
  # define N, should get it from pheno file...
  

  n1=$(grep -e ^\"${p}\" ${pheno_file} | grep -v " ${p}" | awk -F',"' '{print $9}' | sed 's/"//g')
  n2=$(grep -e ^\"${p}\" ${pheno_file} | grep -v " ${p}" | awk -F',"' '{print $10}' | sed 's/"$//g' | sed 's/"//g' | sed 's/$//g')
  
  N1=$(( ${n1} + ${n2} ))
  
  if [ -f ${p_file} ]
   then
   echo 'Run munge_sumstats to create input files for LDSC - ' ${p}
   munge_sumstats.py \
    --sumstats ${p_file} \
    --N ${N1} \
    --p P \
    --snp rsid \
    --out ${p} \
    --a1 a_2 --a2 a_1 \
    --merge-alleles ${ldscores_dir}/eur_w_ld_chr/w_hm3.snplist
  else
   echo Please generate ${p_file} by running prepare_publicGWAS_QC_UKB.sh
  fi
fi


