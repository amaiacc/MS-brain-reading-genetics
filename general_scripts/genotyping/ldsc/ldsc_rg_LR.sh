#!/bin/bash
# PRS
#$ -N ldsc_rg_lat
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
name=${1}
p1=${2}
p2=${3}
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
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/ldsc/UKB_BIG40/
cd ${working_dir}
#--------

if [ ! -f ${working_dir}/${name}_LR_rg.log ]
 then
 if [ -f ${working_dir}/${p1}.sumstats.gz ] && [ -f ${working_dir}/${p2}.sumstats.gz ] 
 then
   echo 'Run LDSC to calculate genetic correlation between left and right for ' ${name}': '${p1}' and '${p2}
   ldsc.py \
   --rg ${p1}.sumstats.gz,${p2}.sumstats.gz, \
   --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --out ${name}_LR_rg
  else
   echo Please generate ${p1}.sumstats.gz and ${p2}.sumstats.gz first.
  fi
  else
  echo File ${name}_LR_rg.log already exists so... rg has probably been computed already.
fi
