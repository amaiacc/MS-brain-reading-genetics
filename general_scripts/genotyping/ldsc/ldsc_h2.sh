#!/bin/bash
# PRS
#$ -N ldsc_h2
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

if [ ! -f ${working_dir}/${p}_h2.log ]
 then
 if [ -f ${working_dir}/${p}.sumstats.gz ]
 then
   echo 'Estimate heritability for' ${p}
   ldsc.py \
   --h2 ${working_dir}/${p}.sumstats.gz \
   --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --out ${working_dir}/${p}_h2
  else
   echo Please generate ${p}.sumstats.gz first.
  fi
  else
  echo File ${p}_h2.log already exists so... h2 has probably been computed already.
fi


