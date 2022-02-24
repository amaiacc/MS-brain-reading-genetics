#!/bin/bash
# Run REML
#$ -N REML_ABCD
#$ -cwd
#$ -q short.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m as

#----------------------------------------------------------------------
# Define paths and directories
#----------------------------------------------------------------------
resource_dir=/export/home/acarrion/acarrion/projects/resources/bin/
gcta_dir=${resource_dir}gcta_1.93.2beta/
PATH=${gcta_dir}:${PATH}
#----------------------------------------------------------------------
# working directory
project_dir=/export/home/acarrion/acarrion/projects/resources/datasets/ABCD/
input_dir=${project_dir}data/working_data/genotyping/preImputation_QC/clean/
working_dir=${project_dir}data/working_data/genotyping/GCTA/
mkdir -p ${working_dir}
mkdir -p ${working_dir}reml/logs
#----------------------------------------------------------------------

#root=ABCD_release_3.0_QCed_EUR.EUR6sd.cleaned # 
#phenotype=nihtbx_reading_uncorrected # 

#----------------------------------------------------------------------
# get parameters from command line arguments:
root=$1
phenotype=$2
pheno_file=$3 # define phenofile here, otherwise definition changes all the time
#pheno_file=${working_dir}/pheno/${phenotype}.txt
#----------------------------------------------------------------------

#----------------------------------------------------------------------
dos2unix $pheno_file
head -n 1 $pheno_file | sed 's/ /\n/g' | sed 's/\t/\n/g' > ${phenotype}.header
# get which column does the phenotype correspond to
pheno_n=$(grep -nr $phenotype'$' ${phenotype}.header | awk -F':' '{print $1}')

# make temporary phenotype file including only : ID1, ID2, phenotype; and no header!
awk -v c1=$pheno_n '{print $1,$2,$c1}' $pheno_file | tail -n +2 | grep -v NA > ${working_dir}/${phenotype}.pheno
echo 'Check header of phenotype file'
head ${working_dir}/${phenotype}.pheno # just to make sure that it's not the same for all

if [ ! -f ${working_dir}/reml/${phenotype}.hsq ]
then
if [ -f ${working_dir}/${phenotype}.pheno ]
then
echo 'One-grm, quantitative, residuals after adjusting for covariates'
gcta64 --reml --grm ${working_dir}/grm/${root}  --pheno ${working_dir}/${phenotype}.pheno --out ${working_dir}/reml/${phenotype}  > ${working_dir}/reml/logs/${phenotype}'.log'
fi
fi
# --thread-num 10 \

# clean intermediate files
rm ${phenotype}.header
rm ${working_dir}/${phenotype}.pheno


