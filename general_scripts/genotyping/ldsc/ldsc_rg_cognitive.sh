#!/bin/bash
# PRS
#$ -N ldsc_rg_cognitive
#$ -cwd
#$ -q short.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m as
#--------
# Run genetic correlation between two traits using LDSC
#--------

#--------
# get arguments
#--------
p1=${1}
gwas_dir=${2}
#--------
name=${p1}_cognitive


#------------------------------------
# GenCor with GWAS summary statistics
#------------------------------------ 
# define phenotypes
p2=Lee2018_EA3
p3=Lee2018_CP
p4=DDyslexia_Gialluisi2020


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
#    --rg ${p1}.sumstats.gz,../${p2}.sumstats.gz,../${p3}.sumstats.gz,../${p4}.sumstats.gz,../${p5}.sumstats.gz,../${p6}.sumstats.gz,../${p7}.sumstats.gz,../${p8}.sumstats.gz,../${p9}.sumstats.gz,../${p10}.sumstats.gz \

if [ ! -f ${working_dir}/${name}_rg.log ]
 then
 if [ -f ${working_dir}/${p1}.sumstats.gz ]
 then
   echo 'Run LDSC to calculate genetic correlation between: '${p1}' and cognitive traits of interest.'
   ldsc.py \
   --rg ${p1}.sumstats.gz,../${p2}.sumstats.gz,../${p3}.sumstats.gz,../${p4}.sumstats.gz \
   --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --out ${name}_rg
  else
   echo Please generate ${p1}.sumstats.gz first.
  fi
  else
  echo File ${name}_rg.log already exists so... rg has probably been computed already.
fi
