#!/bin/bash
# Post-imputation QC per dose file
#$ -N filter_imputed
#$ -cwd
#$ -q long.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m beas
#----------------------------------------------------------------------

module load plink/6.15

#----------------------------------------------------------------------
# Filter imputed genotype vcf files
#----------------------------------------------------------------------
# Once imputation has been ran, we need to filter the SNPs for r2 > 0.8 and MAF > .01. 

#----------------------------------------------------------------------
# Define variables
#----------------------------------------------------------------------
# required parameters (hard-coded or to be taken as arguments)
# define them as command line arguments:
resource_dir=/export/home/acarrion/acarrion/projects/resources/reference_data/
#project_dir=/export/home/acarrion/acarrion/projects/
primary_dir=${1}
working_dir=${2}
data=${3}
i=${4}
#----------------------------------------------------------------------


file=${i}.dose.vcf.gz

#--cd ${working_dir}/tmp--------------------------------------------------------------------
mkdir -p ${working_dir}/tmp ${working_dir}/clean/


if [ -f ${primary_dir}/${file} ]
then

if [ ! -f ${working_dir}/clean/${data}.nodup.r2.bim ] 
then

echo Convert imputed data to plink and filter

if [ -f ${primary_dir}/${i}.info.gz ]; then gunzip ${primary_dir}/${i}.info.gz; fi

if [ ! -f ${i}.bim ]
  then 
    echo "Running ${i}"
    # only calls with genotype probability > 0.9
    plink --double-id --vcf ${primary_dir}${file} --vcf-min-gp 0.90 --biallelic-only --make-bed --out ${i}
    # list of all the SNPs on the chromosome
    plink --bfile ${i} --write-snplist --out ${i}_snps
    # list of duplicates based on the snplist
    cat ${i}_snps.snplist | sort | uniq -d > ${i}_snps.duplicated
    # exclude markers that are duplicated
    plink --bfile ${i} --exclude ${i}_snps.duplicated --make-bed --out ${i}.nodup 
fi

# final bfile set for the chr filtered by r2, based on the info file
plink --bfile ${i}.nodup --qual-scores ${primary_dir}${i}.info 7 1 1 --qual-threshold 0.7 --make-bed --out ${i}.nodup.r2


fi
fi
