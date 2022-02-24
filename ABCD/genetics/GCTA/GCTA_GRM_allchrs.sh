#!/bin/bash
# Build GRM
#$ -N GRM_build
#$ -cwd
#$ -q long.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m beas
#----------------------------------------------------------------------
# get parameters, from command line arguments
root=$1
# root=ABCD_release_3.0_QCed_EUR.EUR6sd.cleaned
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Define paths and directories
#----------------------------------------------------------------------
resource_dir=/export/home/acarrion/acarrion/projects/resources/bin/
gcta_dir=${resource_dir}gcta_1.93.2beta/
PATH=${gcta_dir}:${PATH}
# working directory
project_dir=/export/home/acarrion/acarrion/projects/resources/datasets/ABCD/
input_dir=${project_dir}data/working_data/genotyping/preImputation_QC/clean/
working_dir=${project_dir}data/working_data/genotyping/GCTA/grm/
mkdir -p ${working_dir}

#----------------------------------------------------------------------
mkdir -p ${working_dir} ${working_dir}logs
cd ${working_dir}


echo Generate GRM for all autosomes

# Only run if output file does not exist
if [ ! -f ${root}'.grm.N.bin' ]
then 
 gcta64 --bfile ${input_dir}${root} --maf 0.01 --make-grm  --out ${root} --thread-num 10 > ${working_dir}logs/'grm_'${root}'.log'
fi

#if [ ! -f '/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v2/cal/h2//grm/${root}_rm025_adj.grm.N.bin' ]
#then
#echo 'Prune the GRM by a cutoff of 0.025 (to remove cryptic relatedness): --grm-cutoff  0.025'
#echo 'and adjust for prediction errors assuming the causal variants have similar distribution of allele frequencies as the genotyped SNPs): --grm-adj  0'
#gcta64 --grm ${root} --grm-adj 0 --grm-cutoff  0.025 --make-grm --out ${root}_rm025_adj > ${working_dir}logs/'grm_rm025_adj_'${root}'.log'
#fi


if [ ! -f '/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v2/cal/h2//grm/${root}_rm05_adj.grm.N.bin' ]
then
echo 'Prune the GRM by a cutoff of 0.05 (to remove cryptic relatedness): --grm-cutoff  0.05'
echo 'and adjust for prediction errors assuming the causal variants have similar distribution of allele frequencies as the genotyped SNPs): --grm-adj  0'
gcta64 --grm ${root} --grm-adj 0 --grm-cutoff  0.05 --make-grm --out ${root}_rm05_adj > ${working_dir}logs/'grm_rm05_adj_'${root}'.log'
fi


# use family based data to compute snp h2 and ped h2
#if [ ! -f '/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v2/cal/h2//grm/${root}_bK05_adj.grm.N.bin' ]
#then
#echo 'set the GRM off-diagonal elements that are below the threshold to 0: --make-bK 0.05'
#gcta64 --grm ${root} --make-bK 0.05 --make-grm --out ${root}_bK05_adj > ${working_dir}logs/'grm_bK05_adj_'${root}'.log'
#fi
