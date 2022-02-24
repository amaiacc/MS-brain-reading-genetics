#!/bin/bash
# Run REML bivar
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
phenotype1=$2
phenotype2=$3
pheno_file1=$4 # define phenofile here, otherwise definition changes all the time
pheno_file2=$5
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# get header
head -n 1 $pheno_file1 | sed 's/ /\n/g' | sed 's/\t/\n/g' > ${phenotype1}.header
# get which column does the phenotype correspond to
pheno1_n=$(grep -nr $phenotype1'$' ${phenotype1}.header | awk -F':' '{print $1}')
# make temporary phenotype file including only : ID1, ID2, phenotype; and no header!
awk -v c1=$pheno1_n '{print $1,$2,$c1}' $pheno_file1 | tail -n +2 | grep -v NA > ${working_dir}/${phenotype1}.pheno

# same for pheno2
head -n 1 $pheno_file2 | sed 's/ /\n/g' | sed 's/\t/\n/g' > ${phenotype2}.header
pheno2_n=$(grep -nr $phenotype2'$' ${phenotype2}.header | awk -F':' '{print $1}')
awk -v c1=$pheno2_n '{print $1,$2,$c1}' $pheno_file2 | tail -n +2 | grep -v NA > ${working_dir}/${phenotype2}.pheno

# combine pheno file for both
awk 'FNR==NR{a[$1,$2]=$3;next}{ print $0, a[$1,$2]}' ${working_dir}/${phenotype1}.pheno ${working_dir}/${phenotype2}.pheno | \
 grep -v -e " $" -e "  " > ${working_dir}/${phenotype1}_${phenotype2}.pheno

echo 'Check header of phenotype file'
head ${working_dir}/${phenotype1}_${phenotype2}.pheno # just to make sure that it's not the same for all

if [ ! -f ${working_dir}/reml/${phenotype1}_${phenotype2} ]
then
if [ -f ${working_dir}/${phenotype1}_${phenotype2}.pheno ]
then
echo 'One-grm, quantitative, residuals after adjusting for covariates'
gcta64 --reml --grm ${working_dir}/grm/${root}  --pheno ${working_dir}/${phenotype1}_${phenotype2}.pheno --out ${working_dir}/reml/${phenotype1}_${phenotype2}  > ${working_dir}/reml/logs/${phenotype1}_${phenotype2}'.log'
fi
fi
# --thread-num 10 \

# clean intermediate files
rm ${phenotype1}.header ${phenotype2}.header
rm ${working_dir}/${phenotype1}.pheno  ${working_dir}/${phenotype2}.pheno ${working_dir}/${phenotype1}_${phenotype2}.pheno


